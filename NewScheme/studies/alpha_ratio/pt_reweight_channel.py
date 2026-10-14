"""Independent category x lepton-channel DY hadronic-pT calibrations."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from scipy.stats import chi2 as chi2_distribution
from pt_reweight_core import (Space, Estimate, TAGS, fit_pt, group_edges,
                             merge_bad_groups, predictions, compare)

CALIBRATIONS = {'inclVR': [0, 2], 'heldoutVR': [0]}
CHANNELS = {'ee': [0], 'mumu': [1], 'combined': [0, 1]}


def prefix(cat, tag, channel, target, calibration, axis):
    return '__'.join((cat, tag, channel, target, calibration, axis))+'__'


def add_channels(parts):
    result = {}
    for name in parts[0]:
        if name.startswith('transfer'):
            continue
        result[name] = parts[0][name]
        for part in parts[1:]:
            result[name] = result[name]+part[name]
    return result


def channel_predictions(space, tag, channels, region, groups, fits, cr_regions):
    return add_channels([predictions(space, tag, [channel], region, groups,
                                    fits[channel], cr_regions=cr_regions)
                         for channel in channels])


def channel_mass_groups(space, tag, fits, initial=None):
    """Common mass grid; positivity uses each channel's own factors, DY only."""
    groups = ([np.array([i]) for i in range(space.shape[-2])] if initial is None
              else [g.copy() for g in initial])
    while True:
        bad = np.zeros(len(groups), dtype=bool)
        for channel in (0, 1):
            for region in (0, 1, 2):
                bad |= space.linear(1, [tag], [channel], [region], groups).value <= 0
                for cal in CALIBRATIONS:
                    fit = fits[cal][channel]
                    item = space.weighted([tag], [channel], [region], groups,
                                          fit['groups'], fit['factor'])
                    bad |= item.value <= 0
        if not bad.any():
            return groups
        if len(groups) == 1:
            raise ValueError('Nonpositive integrated DY in a channel/category')
        merge_bad_groups(groups, bad)


def common_pt_groups(fits, channels, npt):
    """Coarsest shared boundaries: never split any channel calibration bin."""
    boundaries = set(range(npt+1))
    for channel in channels:
        boundaries &= set([int(g[0]) for g in fits[channel]['groups']]+[npt])
    edges = sorted(boundaries)
    return [np.arange(a, b) for a, b in zip(edges[:-1], edges[1:])]


def aggregate(item, groups):
    return Estimate(np.array([item.value[g].sum() for g in groups]),
                    np.stack([item.jac[g].sum(axis=0) for g in groups]))


def pt_predictions(space, tag, channels, regions, groups, fits):
    fine = [np.array([i]) for i in range(space.shape[-1])]
    parts = []
    for channel in channels:
        fit = fits[channel]
        obs = space.linear(0, [tag], [channel], regions, ptgroups=groups)
        dy = space.linear(1, [tag], [channel], regions, ptgroups=fine)
        lookup = np.empty(len(fine), dtype=int)
        for i, group in enumerate(fit['groups']):
            lookup[group] = i
        factor = Estimate(fit['factor'].value[lookup], fit['factor'].jac[lookup])
        nominal = aggregate(dy, groups)
        norm = fit['observation'].total()/fit['dy'].total()
        norm = Estimate(np.repeat(norm.value, len(groups)), np.repeat(norm.jac, len(groups), axis=0))
        parts.append(dict(observed=obs, dy_nominal=nominal,
                          dy_reweighted=aggregate(dy*factor, groups), dy_normonly=nominal*norm))
    return add_channels(parts)


def stable_bins(space, tag, channels, regions, groups, axis, nominal):
    var = np.diag(nominal.covariance(space.variance))
    neff = np.divide(nominal.value**2, var, out=np.zeros_like(var), where=var > 0)
    counts = []
    for group in groups:
        masses = group if axis == 'mass' else np.arange(space.shape[-2])
        pts = group if axis == 'pt' else np.arange(space.shape[-1])
        counts.append(space.data[np.ix_([tag], channels, regions, masses, pts)].sum())
    return (np.array(counts) >= 5) & (neff >= 5)


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--min-effective', type=float, default=25.)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError('Use a fresh or empty output directory')
    args.output.mkdir(parents=True, exist_ok=True)
    with np.load(args.source/'source_histograms.npz') as cached:
        arrays = {key: cached[key] for key in cached.files}
    saved, rows, binrows, sfrows, checks, fitrows, projection_checks = {}, [], [], [], [], [], []
    changed_mass_grids = []
    grid_file = args.source/'baseline_mass_edges.npz'
    if not grid_file.exists():
        grid_file = args.source/'results.npz'
    with np.load(grid_file) as archive:
        baseline_edges = {f'{cat}__{tag}__mass_edges': archive[f'{cat}__{tag}__mass_edges']
                          for cat in ('resolved', 'merged') for tag in TAGS}
    previous = {}
    baseline = args.source/'baseline_summary.csv'
    if not baseline.exists():
        baseline = args.source/'summary.csv'
    with baseline.open() as handle:
        for row in csv.DictReader(handle):
            if row['scope'] == 'category':
                previous[tuple(row[k] for k in ('category', 'tag', 'channel', 'target', 'calibration', 'prediction'))] = row

    def store(space, cat, tag, ch, target, cal, axis, edges, items, mask):
        root = prefix(cat, tag, ch, target, cal, axis)
        saved[root+'edges'] = edges
        saved[root+'chi2_mask'] = mask
        observation = items['observed']
        for name, item in items.items():
            cov = item.covariance(space.variance)
            saved[root+name] = np.stack([item.value, np.diag(cov)])
            saved[root+name+'_cov'] = cov
            if name == 'observed':
                continue
            residual = observation-item
            saved[root+name+'_residual_cov'] = residual.covariance(space.variance)
            ratio = np.full((2, len(item.value)), np.nan)
            good = item.value > 0
            if good.any():
                rr = Estimate(observation.value[good], observation.jac[good])/Estimate(item.value[good], item.jac[good])
                ratio[0, good] = rr.value
                ratio[1, good] = np.diag(rr.covariance(space.variance))
            saved[root+name+'_closure_bins'] = ratio
            result = compare(observation, item, space.variance, valid=mask)
            result['chi2_pvalue'] = float(chi2_distribution.sf(result['chi2'], result['rank'])) if result['rank'] else np.nan
            rows.append(dict(category=cat, tag=tag, channel=ch, target=target,
                             calibration=cal, axis=axis, prediction=name, **result))
            for i in range(len(item.value)):
                binrows.append(dict(category=cat, tag=tag, channel=ch, target=target,
                    calibration=cal, axis=axis, prediction=name, low=edges[i], high=edges[i+1],
                    observed=observation.value[i], predicted=item.value[i], variance=cov[i, i],
                    closure=ratio[0, i], closure_variance=ratio[1, i], chi2_used=bool(mask[i])))
            if axis == 'mass' and name in ('dy_nominal', 'alpha_nominal', 'alpha_CRonly_nominal'):
                oldkey = ('alpha_inclusiveCR_nominal' if target == 'SR' and name == 'alpha_nominal'
                          else 'alpha_nominal' if name == 'alpha_CRonly_nominal' else name)
                old = previous[(cat, tag, ch, target, cal, oldkey)]
                fields = ('observed', 'observed_stat', 'predicted', 'predicted_stat')
                if (cat, tag) in changed_mass_grids and name != 'dy_nominal':
                    fields = ('observed', 'observed_stat')
                for field in fields:
                    np.testing.assert_allclose(result[field], float(old[field]), rtol=1e-10, atol=1e-7)
        return items

    for cat in ('resolved', 'merged'):
        space = Space(arrays, cat)
        for ti, tag in enumerate(TAGS):
            fits = {cal: {channel: fit_pt(space, [ti], regions, channels=[channel],
                                         min_effective=args.min_effective)
                          for channel in (0, 1)} for cal, regions in CALIBRATIONS.items()}
            for cal, bychannel in fits.items():
                for channel, fit in bychannel.items():
                    ch = 'ee' if channel == 0 else 'mumu'
                    root = f'{cat}__{tag}__{ch}__SF__{cal}__'
                    covariance = fit['factor'].covariance(space.variance)
                    saved[root+'edges'] = fit['edges']
                    saved[root+'values'] = np.stack([fit['factor'].value, np.diag(covariance)])
                    saved[root+'cov'] = covariance
                    fitrows.append(dict(category=cat, tag=tag, channel=ch, calibration=cal,
                                        pt_bins=len(fit['groups']), pt_edges=fit['edges'].tolist()))
                    other_indices = np.r_[space.indices(0, [ti], [1-channel], [0, 1, 2]),
                                          space.indices(1, [ti], [1-channel], [0, 1, 2])]
                    assert np.all(fit['factor'].jac[:, other_indices] == 0)
                    assert np.all(fit['factor'].jac[:, space.indices(0, [ti], [channel], [1])] == 0)
                    for i in range(len(fit['groups'])):
                        sfrows.append(dict(category=cat, tag=tag, channel=ch, calibration=cal,
                            pt_low=fit['edges'][i], pt_high=fit['edges'][i+1], sf=fit['factor'].value[i],
                            sf_stat=np.sqrt(covariance[i, i]), effective_Z=fit['effective_z'][i],
                            effective_DY=fit['effective_dy'][i]))
            old_edges = baseline_edges[f'{cat}__{tag}__mass_edges']
            old_indices = np.searchsorted(space.mass_edges, old_edges)
            np.testing.assert_array_equal(space.mass_edges[old_indices], old_edges)
            initial = [np.arange(a, b) for a, b in zip(old_indices[:-1], old_indices[1:])]
            mg = channel_mass_groups(space, ti, fits, initial=initial)
            medges = group_edges(space.mass_edges, mg)
            if not np.array_equal(medges, old_edges):
                changed_mass_grids.append((cat, tag))
            for cal, regions in CALIBRATIONS.items():
                for ch, channels in CHANNELS.items():
                    masses = {}
                    for target, region in [('SR', 1), ('VR', 2)]:
                        result = channel_predictions(space, ti, channels, region, mg, fits[cal],
                                                     cr_regions=(0, 2) if target == 'SR' else (0,))
                        if target == 'SR':
                            alt = channel_predictions(space, ti, channels, region, mg, fits[cal], cr_regions=(0,))
                            result['alpha_CRonly_nominal'] = alt['alpha_nominal']
                            result['alpha_CRonly_reweighted'] = alt['alpha_reweighted']
                        mask = stable_bins(space, ti, channels, [region], mg, 'mass', result['dy_nominal'])
                        masses[target] = store(space, cat, tag, ch, target, cal, 'mass', medges, result, mask)
                    pg = common_pt_groups(fits[cal], channels, space.shape[-1])
                    pedges = group_edges(space.pt_edges, pg)
                    for target, rr in [('SR', [1]), ('calibration', regions)]:
                        result = pt_predictions(space, ti, channels, rr, pg, fits[cal])
                        mask = stable_bins(space, ti, channels, rr, pg, 'pt', result['dy_nominal'])
                        store(space, cat, tag, ch, target, cal, 'pt', pedges, result, mask)
                        if target == 'calibration':
                            np.testing.assert_allclose(result['dy_reweighted'].value, result['observed'].value,
                                                       rtol=1e-12, atol=1e-9)
                            np.testing.assert_allclose(result['dy_reweighted'].jac, result['observed'].jac,
                                                       rtol=0, atol=1e-12)
                            checks.append(dict(category=cat, tag=tag, channel=ch, calibration=cal,
                                max_absolute_residual=float(np.max(np.abs(result['observed'].value-result['dy_reweighted'].value)))))
                        else:
                            for name in result:
                                left, right = result[name].total(), masses['SR'][name].total()
                                np.testing.assert_allclose(left.value, right.value, rtol=1e-12, atol=1e-8)
                                np.testing.assert_allclose(left.jac, right.jac, rtol=0, atol=1e-12)
                            projection_checks.append(f'{cat}/{tag}/{ch}/{cal}')
            print(f'{cat} {tag}: independent ee/mumu fits, {len(mg)} mass bins; closure and projections passed', flush=True)
    write_csv(args.output/'summary.csv', rows)
    write_csv(args.output/'bins.csv', binrows)
    write_csv(args.output/'pt_factors.csv', sfrows)
    np.savez_compressed(args.output/'results.npz', **saved)
    metadata = dict(date='2026-10-07', source=str(args.source.resolve()),
        source_sha256=hashlib.sha256((args.source/'source_histograms.npz').read_bytes()).hexdigest(),
        correction='Independent absolute (Data-TTplusWW-VZ)/DY in each topology/tag/lepton channel',
        shared_tag_corrections=False, ZH_subtracted=False, min_effective_pt=args.min_effective,
        combined_channels='Sum predictions from separate ee/mumu corrections and transfers; no combined fitted correction',
        calibrations=CALIBRATIONS, alpha_source='SR:C+VR; VR:C disjoint from target; SR C-only sensitivity also saved',
        statistics='Full first-order joint covariance, calibration reuse retained; separate-channel fits use disjoint primitives',
        missing_systematics='Lumi, tag SFs, subtraction/modeling, MC normalization covariance and adaptive-bin-choice uncertainty',
        chi2='Stat-only Gaussian full covariance on data count>=5 and DY Neff>=5; not final likelihood',
        fits=fitrows, calibration_checks=checks, mass_pt_projection_checks=projection_checks,
        nominal_yield_and_error_comparison='Same baseline mass grid unless additional DY-only merging required; unchanged-grid nominal predictions checked',
        additionally_merged_mass_grids=changed_mass_grids,
        mass_bins='DY-only positivity of each channel using its own calibration; no target observations in bin choice',
        combined_pt_bins='Common coarsening of channel fit grids; each channel factor remains unchanged')
    (args.output/'metadata.json').write_text(json.dumps(metadata, indent=2))
    for name in ('source_histograms.npz', 'source_audit.json'):
        shutil.copy2(args.source/name, args.output/name)
    shutil.copy2(baseline, args.output/'baseline_summary.csv')
    np.savez_compressed(args.output/'baseline_mass_edges.npz', **baseline_edges)
    print('Numerical channel-dependent analysis complete: '+str(args.output), flush=True)


if __name__ == '__main__':
    main()
