"""Calibration-region mass shapes before/after category/channel pT correction."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import zipfile
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import mplhep as hep
import numpy as np
from scipy.stats import chi2 as chi2_distribution
from pt_reweight_core import Space, Estimate, TAGS, fit_pt, compare
from pt_reweight_channel import (CALIBRATIONS, CHANNELS, prefix, add_channels,
                                stable_bins, write_csv)
from plot_pt_reweight_channel import comparison


def direct_mass(space, tag, channels, regions, groups, fits):
    parts = []
    for channel in channels:
        fit = fits[channel]
        obs = space.linear(0, [tag], [channel], regions, groups)
        dy = space.linear(1, [tag], [channel], regions, groups)
        rw = space.weighted([tag], [channel], regions, groups, fit['groups'], fit['factor'])
        norm = fit['observation'].total()/fit['dy'].total()
        norm = Estimate(np.repeat(norm.value, len(groups)), np.repeat(norm.jac, len(groups), axis=0))
        parts.append(dict(observed=obs, dy_nominal=dy, dy_reweighted=rw, dy_normonly=dy*norm))
    return add_channels(parts)


def figure(data, summaries, cat, tag, ch, target, cal, compact=False, axes=None):
    fig = comparison(data, summaries, cat, tag, ch, target, cal, 'mass', 'direct', compact=compact, axes=axes)
    root = prefix(cat, tag, ch, target, cal, 'mass')
    lower = fig.axes[-1] if axes is None else axes[1]
    from plot_pt_reweight import display_edges
    edges = display_edges(data[root+'edges'])
    center = (edges[:-1]+edges[1:])/2
    ratio, var = data[root+'dy_normonly_closure_bins']
    mask = data[root+'chi2_mask'].astype(bool) & np.isfinite(ratio)
    lower.errorbar(center[mask], ratio[mask], yerr=np.sqrt(np.maximum(0., var[mask])),
                   fmt='^', ms=2 if compact else 2.5, color='#666666', elinewidth=.7)
    lo, hi = lower.get_ylim()
    if mask.any():
        lower.set_ylim(min(lo, np.min(ratio[mask])-.1), max(hi, np.max(ratio[mask])+.1))
    if not compact:
        for item in fig.texts:
            if abs(item.get_position()[1]-.055) < .001:
                item.set_text('In-sample check: pT closure enforced; mass-shape closure not enforced.')
        metrics = []
        for key, label in [('dy_nominal', 'nominal'), ('dy_normonly', 'norm-only'), ('dy_reweighted', 'pT')]:
            row = summaries[(cat, tag, ch, target, cal, 'mass', key)]
            metrics.append(f"{label} {row['chi2']:.1f}/{row['rank']}")
        fig.text(.5, .12, 'Mass '+r'$\chi^2$/rank: '+', '.join(metrics), ha='center', fontsize=10)
    return fig


def refresh_package(folder):
    assert len(list(folder.glob('*.png'))) == 42
    assert len(list(folder.glob('*.pdf'))) == 43
    assert all(p.with_suffix('.pdf').exists() for p in folder.glob('*.png'))
    validation = json.loads((folder/'validation.json').read_text())
    assert validation['status'] == 'passed'
    digest = hashlib.sha256((folder/'source_histograms.npz').read_bytes()).hexdigest()
    assert digest == validation['source_sha256']
    doc = Path(__file__).parent/'CR_MASS_PT_CHANNEL.md'
    if doc.exists():
        shutil.copy2(doc, folder/'README.md')
    if Path(__file__).resolve().parent != folder.resolve():
        shutil.copy2(Path(__file__), folder/Path(__file__).name)
    validation['script_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in folder.glob('*.py')}
    (folder/'validation.json').write_text(json.dumps(validation, indent=2))
    with zipfile.ZipFile(folder.with_suffix('.zip'), 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(folder.iterdir()):
            if path.is_file():
                archive.write(path, arcname=folder.name+'/'+path.name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--finalize-only', action='store_true')
    args = parser.parse_args()
    if args.finalize_only:
        refresh_package(args.output)
        print('Verified 42 PNG/PDF pairs and refreshed documentation/package.', flush=True)
        return
    args.output.mkdir(parents=True, exist_ok=False)
    reference_path = args.study/'calibration_reference.npz'
    if not reference_path.exists():
        reference_path = args.study/'results.npz'
    metadata_path = args.study/'reference_metadata.json'
    if not metadata_path.exists():
        metadata_path = args.study/'metadata.json'
    reference_metadata = json.loads(metadata_path.read_text())
    digest = hashlib.sha256((args.study/'source_histograms.npz').read_bytes()).hexdigest()
    assert digest == reference_metadata['source_sha256']
    with np.load(args.study/'source_histograms.npz') as cached:
        arrays = {key: cached[key] for key in cached.files}
    with np.load(reference_path) as archive:
        reference = {key: archive[key] for key in archive.files
                     if '__SF__' in key or key.endswith('__edges')}
    saved, rows, bins, checks, summaries = {}, [], [], [], {}
    for cat in ('resolved', 'merged'):
        space = Space(arrays, cat)
        for ti, tag in enumerate(TAGS):
            edges = reference[prefix(cat, tag, 'ee', 'SR', 'inclVR', 'mass')+'edges']
            indices = np.searchsorted(space.mass_edges, edges)
            np.testing.assert_array_equal(space.mass_edges[indices], edges)
            groups = [np.arange(a, b) for a, b in zip(indices[:-1], indices[1:])]
            for cal, regions in CALIBRATIONS.items():
                target = 'CR+VR' if cal == 'inclVR' else 'CR'
                fits = {channel: fit_pt(space, [ti], regions, channels=[channel],
                                       min_effective=reference_metadata['min_effective_pt'])
                        for channel in (0, 1)}
                for channel, fit in fits.items():
                    ch = 'ee' if channel == 0 else 'mumu'
                    sfroot = f'{cat}__{tag}__{ch}__SF__{cal}__'
                    np.testing.assert_array_equal(fit['edges'], reference[sfroot+'edges'])
                    np.testing.assert_allclose(fit['factor'].value, reference[sfroot+'values'][0], rtol=1e-12)
                    np.testing.assert_allclose(fit['factor'].covariance(space.variance), reference[sfroot+'cov'], rtol=1e-10, atol=1e-12)
                for ch, channels in CHANNELS.items():
                    items = direct_mass(space, ti, channels, regions, groups, fits)
                    obs = items['observed']
                    for key in ('dy_reweighted', 'dy_normonly'):
                        np.testing.assert_allclose(items[key].total().value, obs.total().value, rtol=1e-12, atol=1e-8)
                        np.testing.assert_allclose(items[key].total().jac, obs.total().jac, rtol=0, atol=1e-12)
                    checks.append(f'{cat}/{tag}/{ch}/{cal}')
                    mask = stable_bins(space, ti, channels, regions, groups, 'mass', items['dy_nominal'])
                    root = prefix(cat, tag, ch, target, cal, 'mass')
                    saved[root+'edges'], saved[root+'chi2_mask'] = edges, mask
                    for key, item in items.items():
                        covariance = item.covariance(space.variance)
                        saved[root+key] = np.stack([item.value, np.diag(covariance)])
                        saved[root+key+'_cov'] = covariance
                        if key == 'observed':
                            continue
                        result = compare(obs, item, space.variance, valid=mask)
                        result['chi2_pvalue'] = float(chi2_distribution.sf(result['chi2'], result['rank'])) if result['rank'] else np.nan
                        row = dict(category=cat, tag=tag, channel=ch, target=target, calibration=cal, axis='mass', prediction=key, **result)
                        rows.append(row)
                        summaries[(cat, tag, ch, target, cal, 'mass', key)] = row
                        saved[root+key+'_residual_cov'] = (obs-item).covariance(space.variance)
                        ratio = np.full((2, len(groups)), np.nan)
                        good = item.value > 0
                        if good.any():
                            rr = Estimate(obs.value[good], obs.jac[good])/Estimate(item.value[good], item.jac[good])
                            ratio[0, good], ratio[1, good] = rr.value, np.diag(rr.covariance(space.variance))
                        saved[root+key+'_closure_bins'] = ratio
                        for i in range(len(groups)):
                            bins.append(dict(category=cat, tag=tag, channel=ch, calibration=cal, prediction=key,
                                low=edges[i], high=edges[i+1], observed=obs.value[i], predicted=item.value[i],
                                variance=covariance[i, i], closure=ratio[0, i], closure_variance=ratio[1, i], chi2_used=bool(mask[i])))
    write_csv(args.output/'summary.csv', rows)
    write_csv(args.output/'bins.csv', bins)
    np.savez_compressed(args.output/'results.npz', **saved)
    validation = dict(status='numerical_passed', source_sha256=digest, calibration_factors_unchanged=24,
                      integrated_value_and_jacobian_checks=checks, mass_grid='Same as SR comparison',
                      interpretation='Calibration sample: integrated yield closure enforced, mass shape not enforced; statistical-only in-sample diagnostic')
    (args.output/'validation.json').write_text(json.dumps(validation, indent=2))
    print('Numerical checks passed; rendering calibration mass distributions.', flush=True)
    plt.style.use(hep.style.CMS)
    plt.rcParams.update({'font.size': 14, 'axes.labelpad': 7, 'pdf.fonttype': 42, 'axes.linewidth': 1.1})
    count = 0
    with PdfPages(args.output/'all_CR_mass_pt_reweight.pdf') as bundle:
        def save(fig, stem):
            nonlocal count
            fig.savefig(args.output/(stem+'.pdf'), bbox_inches='tight')
            fig.savefig(args.output/(stem+'.png'), dpi=180, bbox_inches='tight')
            bundle.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            count += 1
        for cal in CALIBRATIONS:
            target = 'CR+VR' if cal == 'inclVR' else 'CR'
            for ch in CHANNELS:
                fig = plt.figure(figsize=(17, 11.5))
                outer = fig.add_gridspec(2, 3, wspace=.32, hspace=.38)
                for i, (cat, tag) in enumerate((c, t) for c in ('resolved', 'merged') for t in TAGS):
                    grid = outer[i//3, i%3].subgridspec(2, 1, height_ratios=[3, 1], hspace=.05)
                    ax = fig.add_subplot(grid[0])
                    lower = fig.add_subplot(grid[1], sharex=ax)
                    figure(saved, summaries, cat, tag, ch, target, cal, compact=True, axes=(ax, lower))
                handles, labels = fig.axes[0].get_legend_handles_labels()
                fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .995), ncol=4, fontsize=11, frameon=False)
                fig.subplots_adjust(top=.90, bottom=.10, left=.075, right=.985)
                fig.text(.5, .015, f'{target} calibration sample; pT closure enforced, mass-shape closure not enforced. Statistical errors only.', ha='center', fontsize=11)
                save(fig, f'{ch}_{cal}_mass_summary')
            for cat in ('resolved', 'merged'):
                for tag in TAGS:
                    for ch in CHANNELS:
                        save(figure(saved, summaries, cat, tag, ch, target, cal), f'{cat}_{tag}_{ch}_{cal}_mass')
    assert count == 42 and len(list(args.output.glob('*.png'))) == 42
    validation.update(status='passed', individual_png_pdf_pairs=count, PDF_files=43, bundle_pages=count)
    (args.output/'validation.json').write_text(json.dumps(validation, indent=2))
    for name in ('source_histograms.npz', 'source_audit.json'):
        shutil.copy2(args.study/name, args.output/name)
    np.savez_compressed(args.output/'calibration_reference.npz', **reference)
    shutil.copy2(metadata_path, args.output/'reference_metadata.json')
    for name in ('plot_cr_mass_pt_reweight_channel.py', 'pt_reweight_channel.py',
                 'plot_pt_reweight_channel.py', 'pt_reweight_core.py', 'plot_pt_reweight.py'):
        shutil.copy2(Path(__file__).parent/name, args.output/name)
    doc = Path(__file__).parent/'CR_MASS_PT_CHANNEL.md'
    if not doc.exists():
        doc = Path(__file__).parent/'README.md'
    shutil.copy2(doc, args.output/'README.md')
    refresh_package(args.output)
    print('Saved 42 CR mass PNG/PDF pairs and 42-page bundle.', flush=True)


if __name__ == '__main__':
    main()
