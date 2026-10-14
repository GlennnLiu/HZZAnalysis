"""Category/channel CR mass reweighting of nominal DY, with joint covariance."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.stats import chi2 as chi2_distribution
from pt_reweight_core import Space, Estimate, TAGS, group_edges, merge_bad_groups, compare
from pt_reweight_channel import CALIBRATIONS, CHANNELS, prefix, add_channels, stable_bins, write_csv


def mass_grid(space, tag, initial, min_effective=25.):
    """Common grid; independent fits. Calibration statistics/target DY, never SR data."""
    groups = [group.copy() for group in initial]
    while True:
        bad = np.zeros(len(groups), dtype=bool)
        for ch in (0, 1):
            for regions in CALIBRATIONS.values():
                z = space.linear(0, [tag], [ch], regions, groups)
                dy = space.linear(1, [tag], [ch], regions, groups)
                zv, dv = np.diag(z.covariance(space.variance)), np.diag(dy.covariance(space.variance))
                nz = np.divide(z.value**2, zv, out=np.zeros_like(zv), where=zv > 0)
                nd = np.divide(dy.value**2, dv, out=np.zeros_like(dv), where=dv > 0)
                bad |= (z.value <= 0) | (dy.value <= 0) | (nz < min_effective) | (nd < min_effective)
            for region in (1, 2):
                bad |= space.linear(1, [tag], [ch], [region], groups).value <= 0
        if not bad.any():
            return groups
        if len(groups) == 1:
            raise ValueError('Insufficient integrated calibration/DY statistics')
        merge_bad_groups(groups, bad)


def fit_mass(space, tag, channel, regions, groups):
    z = space.linear(0, [tag], [channel], regions, groups)
    dy = space.linear(1, [tag], [channel], regions, groups)
    if np.any(z.value <= 0) or np.any(dy.value <= 0):
        raise ValueError('Nonpositive mass calibration')
    return dict(factor=z/dy, observation=z, dy=dy)


def expand_fit(fit, correction_groups, alpha_groups):
    lookup = []
    for group in alpha_groups:
        matches = [i for i, coarse in enumerate(correction_groups) if np.isin(group, coarse).all()]
        if len(matches) != 1:
            raise ValueError('Alpha bin crosses a correction boundary')
        lookup.append(matches[0])
    return dict(fit, factor=Estimate(fit['factor'].value[lookup], fit['factor'].jac[lookup]))


def mass_predictions(space, tag, channels, regions, source_regions, groups, fits):
    if set(regions) & set(source_regions):
        raise ValueError('Alpha source overlaps target')
    parts = []
    for ch in channels:
        fit = fits[ch]
        zt = space.linear(0, [tag], [ch], regions, groups)
        dt = space.linear(1, [tag], [ch], regions, groups)
        zc = space.linear(0, [tag], [ch], source_regions, groups)
        dc = space.linear(1, [tag], [ch], source_regions, groups)
        wt, wc = dt*fit['factor'], dc*fit['factor']
        norm = fit['observation'].total()/fit['dy'].total()
        norm = Estimate(np.repeat(norm.value, len(groups)), np.repeat(norm.jac, len(groups), axis=0))
        parts.append(dict(observed=zt, dy_nominal=dt, dy_reweighted=wt, dy_normonly=dt*norm,
                          alpha_nominal=(dt/dc)*zc, alpha_reweighted=(wt/wc)*zc))
    return add_channels(parts)


def calibration_predictions(space, tag, channels, regions, groups, fits):
    parts = []
    for ch in channels:
        fit = fits[ch]
        norm = fit['observation'].total()/fit['dy'].total()
        norm = Estimate(np.repeat(norm.value, len(groups)), np.repeat(norm.jac, len(groups), axis=0))
        obs = space.linear(0, [tag], [ch], regions, groups)
        dy = space.linear(1, [tag], [ch], regions, groups)
        parts.append(dict(observed=obs, dy_nominal=dy,
                          dy_reweighted=dy*fit['factor'], dy_normonly=dy*norm))
    return add_channels(parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--min-effective', type=float, default=25.)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError('Use a fresh or empty output path')
    args.output.mkdir(parents=True, exist_ok=True)
    with np.load(args.source/'source_histograms.npz') as archive:
        arrays = {key: archive[key] for key in archive.files}
    edge_path = args.source/'baseline_mass_edges.npz'
    if not edge_path.exists():
        edge_path = args.source/'results.npz'
    with np.load(edge_path) as archive:
        baseline = {f'{cat}__{tag}__mass_edges': archive[f'{cat}__{tag}__mass_edges']
                    for cat in ('resolved', 'merged') for tag in TAGS}
    reference_path = args.source/'reference_summary.csv'
    if not reference_path.exists():
        reference_path = args.source/'summary.csv'
    import csv
    with reference_path.open() as handle:
        reference_rows = {tuple(r[k] for k in ('category', 'tag', 'channel', 'target', 'calibration', 'prediction')): r
                          for r in csv.DictReader(handle) if r.get('axis') == 'mass'}
    saved, rows, bins, factors, grids, checks = {}, [], [], [], [], []

    def store(space, cat, tag, ch, target, cal, edges, items, mask):
        root = prefix(cat, tag, ch, target, cal, 'mass')
        saved[root+'edges'], saved[root+'chi2_mask'] = edges, mask
        obs = items['observed']
        for name, item in items.items():
            cov = item.covariance(space.variance)
            saved[root+name] = np.stack([item.value, np.diag(cov)])
            saved[root+name+'_cov'] = cov
            if name == 'observed':
                continue
            saved[root+name+'_residual_cov'] = (obs-item).covariance(space.variance)
            ratio = np.full((2, len(item.value)), np.nan)
            good = item.value > 0
            if good.any():
                rr = Estimate(obs.value[good], obs.jac[good])/Estimate(item.value[good], item.jac[good])
                ratio[0, good], ratio[1, good] = rr.value, np.diag(rr.covariance(space.variance))
            saved[root+name+'_closure_bins'] = ratio
            result = compare(obs, item, space.variance, valid=mask)
            if target in ('SR', 'VR') and name in ('dy_nominal', 'alpha_nominal', 'alpha_CRonly_nominal'):
                old = reference_rows[(cat, tag, ch, target, cal, name)]
                for field in ('observed', 'observed_stat', 'predicted', 'predicted_stat'):
                    np.testing.assert_allclose(result[field], float(old[field]), rtol=1e-10, atol=1e-7)
            result['chi2_pvalue'] = float(chi2_distribution.sf(result['chi2'], result['rank'])) if result['rank'] else np.nan
            rows.append(dict(category=cat, tag=tag, channel=ch, target=target,
                             calibration=cal, axis='mass', prediction=name, **result))
            for i in range(len(item.value)):
                bins.append(dict(category=cat, tag=tag, channel=ch, target=target,
                    calibration=cal, prediction=name, low=edges[i], high=edges[i+1],
                    observed=obs.value[i], predicted=item.value[i], variance=cov[i, i],
                    closure=ratio[0, i], closure_variance=ratio[1, i], chi2_used=bool(mask[i])))

    for cat in ('resolved', 'merged'):
        space = Space(arrays, cat)
        for ti, tag in enumerate(TAGS):
            old_edges = baseline[f'{cat}__{tag}__mass_edges']
            indices = np.searchsorted(space.mass_edges, old_edges)
            np.testing.assert_array_equal(space.mass_edges[indices], old_edges)
            initial = [np.arange(a, b) for a, b in zip(indices[:-1], indices[1:])]
            groups = mass_grid(space, ti, initial, args.min_effective)
            edges = group_edges(space.mass_edges, groups)
            grids.append(dict(category=cat, tag=tag, baseline_edges=old_edges.tolist(), edges=edges.tolist(),
                              alpha_edges=old_edges.tolist(),
                              additionally_merged=not np.array_equal(edges, old_edges)))
            fits = {cal: {ch: fit_mass(space, ti, ch, regions, groups) for ch in (0, 1)}
                    for cal, regions in CALIBRATIONS.items()}
            for cal, regions in CALIBRATIONS.items():
                for ci, fit in fits[cal].items():
                    ch = 'ee' if ci == 0 else 'mumu'
                    root = f'{cat}__{tag}__{ch}__SF__{cal}__'
                    cov = fit['factor'].covariance(space.variance)
                    saved[root+'edges'] = edges
                    saved[root+'values'] = np.stack([fit['factor'].value, np.diag(cov)])
                    saved[root+'cov'] = cov
                    for proc in (0, 1):
                        assert np.all(fit['factor'].jac[:, space.indices(proc, [ti], [1-ci], [0, 1, 2])] == 0)
                    assert np.all(fit['factor'].jac[:, space.indices(0, [ti], [ci], [1])] == 0)
                    if cal == 'heldoutVR':
                        assert np.all(fit['factor'].jac[:, space.indices(0, [ti], [ci], [2])] == 0)
                    for i in range(len(groups)):
                        factors.append(dict(category=cat, tag=tag, channel=ch, calibration=cal,
                            mass_low=edges[i], mass_high=edges[i+1], factor=fit['factor'].value[i],
                            factor_stat=np.sqrt(cov[i, i])))
                for ch, channels in CHANNELS.items():
                    calitems = calibration_predictions(space, ti, channels, regions, groups, fits[cal])
                    np.testing.assert_allclose(calitems['dy_reweighted'].value, calitems['observed'].value, rtol=1e-12, atol=1e-8)
                    np.testing.assert_allclose(calitems['dy_reweighted'].jac, calitems['observed'].jac, rtol=0, atol=1e-12)
                    finefits = {ci: expand_fit(fit, groups, initial) for ci, fit in fits[cal].items()}
                    finecal = calibration_predictions(space, ti, channels, regions, initial, finefits)
                    np.testing.assert_allclose(finecal['dy_reweighted'].total().value, finecal['observed'].total().value, rtol=1e-12, atol=1e-8)
                    np.testing.assert_allclose(finecal['dy_reweighted'].total().jac, finecal['observed'].total().jac, rtol=0, atol=1e-12)
                    mask = stable_bins(space, ti, channels, regions, initial, 'mass', finecal['dy_nominal'])
                    store(space, cat, tag, ch, 'calibration', cal, old_edges, finecal, mask)
                    checks.append(dict(test='mass_calibration_closure', category=cat, tag=tag, channel=ch, calibration=cal))
                    for target, region in [('SR', 1), ('VR', 2)]:
                        source = [0, 2] if target == 'SR' else [0]
                        items = mass_predictions(space, ti, channels, [region], source, initial, finefits)
                        np.testing.assert_allclose(items['alpha_reweighted'].value, items['alpha_nominal'].value, rtol=1e-12, atol=1e-8)
                        np.testing.assert_allclose(items['alpha_reweighted'].jac, items['alpha_nominal'].jac, rtol=0, atol=2e-12)
                        checks.append(dict(test='alpha_factor_cancellation', category=cat, tag=tag, channel=ch, calibration=cal, target=target))
                        if source == regions:
                            coarse = mass_predictions(space, ti, channels, [region], source, groups, fits[cal])
                            np.testing.assert_allclose(coarse['dy_reweighted'].value, coarse['alpha_nominal'].value, rtol=1e-12, atol=1e-8)
                            np.testing.assert_allclose(coarse['dy_reweighted'].jac, coarse['alpha_nominal'].jac, rtol=0, atol=2e-12)
                            checks.append(dict(test='direct_equals_nominal_alpha', category=cat, tag=tag, channel=ch, calibration=cal, target=target))
                        if target == 'SR':
                            alt = mass_predictions(space, ti, channels, [region], [0], initial, finefits)
                            items['alpha_CRonly_nominal'], items['alpha_CRonly_reweighted'] = alt['alpha_nominal'], alt['alpha_reweighted']
                        mask = stable_bins(space, ti, channels, [region], initial, 'mass', items['dy_nominal'])
                        store(space, cat, tag, ch, target, cal, old_edges, items, mask)
            print(f'{cat} {tag}: {len(initial)} -> {len(groups)} mass bins; exact calibration and alpha-cancellation checks passed', flush=True)
    write_csv(args.output/'summary.csv', rows)
    write_csv(args.output/'bins.csv', bins)
    write_csv(args.output/'mass_factors.csv', factors)
    np.savez_compressed(args.output/'results.npz', **saved)
    metadata = dict(date='2026-10-07', correction='Independent topology/tag/lepton-channel binned mass correction from nominal DY; no pT correction',
        source_sha256=hashlib.sha256((args.source/'source_histograms.npz').read_bytes()).hexdigest(),
        calibrations=CALIBRATIONS, ZH_subtracted=False, category_inclusive=False,
        alpha_source='SR:C+VR; VR:C; SR C-only sensitivity saved',
        mass_grid='Start from prior grid; merge until each calibration/channel has Z,DY Neff>=25 and target DY positive; no SR data used; primary calibration includes VR',
        min_effective=args.min_effective, grids=grids, checks=checks,
        cancellation='Correction is constant in every alpha bin; weighted numerator/denominator cancel in values and Jacobians',
        alpha_grid='Original comparison/alpha grid preserved. Coarse mass correction expanded as a constant over its constituent alpha bins; nominal yields/errors checked against earlier study.',
        direct_identity='If mass calibration uses the alpha source, direct weighted DY equals nominal alpha on the correction grid; smoothing within coarser correction bins can differ from fine-bin alpha.',
        uncertainties='Stat-only joint first-order covariance; same-calibration reuse retained',
        limitations='Correction grids can be coarser; plotted mass/alpha grid remains unchanged. Missing tag/lumi/subtraction/modeling/MC-normalization covariance and adaptive-bin-choice uncertainty.')
    (args.output/'metadata.json').write_text(json.dumps(metadata, indent=2))
    for name in ('source_histograms.npz', 'source_audit.json'):
        shutil.copy2(args.source/name, args.output/name)
    np.savez_compressed(args.output/'baseline_mass_edges.npz', **baseline)
    shutil.copy2(reference_path, args.output/'reference_summary.csv')
    print('Mass-reweighting numerical analysis complete.', flush=True)


if __name__ == '__main__':
    main()
