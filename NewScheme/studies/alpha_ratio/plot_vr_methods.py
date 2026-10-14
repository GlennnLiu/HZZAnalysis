"""Validate and compare five VR predictions, with primary and CR-only corrections."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import zipfile
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import mplhep as hep
import numpy as np
from scipy.stats import chi2 as chi2_distribution

TAGS = ('btag', 'untag', 'vbftag')
CHANNELS = ('ee', 'mumu', 'combined')
CALIBRATIONS = ('inclVR', 'heldoutVR')
LABELS = {'btag': 'b-tag', 'untag': 'untagged', 'vbftag': 'VBF-tag'}
CH = {'ee': r'$ee$', 'mumu': r'$\mu\mu$', 'combined': r'$ee+\mu\mu$'}
METHODS = (
    ('direct_dy', 'pt', 'dy_nominal', 'Direct DY', '#e68a00', '-'),
    ('direct_pt', 'pt', 'dy_reweighted', 'Direct DY, pT reweighted', '#2166ac', '--'),
    ('direct_mass', 'mass', 'dy_reweighted', 'Direct DY, mass reweighted', '#b2182b', ':'),
    ('alpha', 'pt', 'alpha_nominal', 'Alpha estimate', '#009E73', '-.'),
    ('alpha_pt', 'pt', 'alpha_reweighted', 'Alpha estimate, pT reweighted', '#8b4daf', (0, (5, 2, 1, 2))),
)


def load(folder):
    with np.load(folder/'results.npz') as archive:
        data = {key: archive[key] for key in archive.files}
    with (folder/'summary.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    summaries = {tuple(row[k] for k in ('category', 'tag', 'channel', 'target',
                 'calibration', 'axis', 'prediction')): row for row in rows}
    return data, summaries, json.loads((folder/'metadata.json').read_text())


def prefix(cat, tag, ch, cal):
    return '__'.join((cat, tag, ch, 'VR', cal, 'mass'))+'__'


def recompute_statistics(observed, predicted, covariance, mask):
    """Same full residual-covariance statistic as the source, no fitted VR norm."""
    matrix = covariance[np.ix_(mask, mask)]
    residual = (observed-predicted)[mask]
    if not matrix.size:
        return np.nan, 0
    eigen, rotation = np.linalg.eigh(matrix)
    keep = eigen > max(float(eigen.max())*1e-10, 1e-12)
    projected = rotation[:, keep].T@residual
    return float(np.sum(projected**2/eigen[keep])), int(keep.sum())


def prepare(pt, mass):
    """Fail closed if nominal inputs, bins or covariance differ across studies."""
    if pt[2]['source_sha256'] != mass[2]['source_sha256']:
        raise ValueError('Different primitive histogram inputs')
    saved, rows = {}, []
    for cal in CALIBRATIONS:
      for cat in ('resolved', 'merged'):
        for tag in TAGS:
            for ch in CHANNELS:
                root = prefix(cat, tag, ch, cal)
                for key in ('edges', 'observed', 'observed_cov', 'chi2_mask',
                            'dy_nominal', 'dy_nominal_cov', 'alpha_nominal', 'alpha_nominal_cov'):
                    np.testing.assert_allclose(pt[0][root+key], mass[0][root+key], rtol=1e-10, atol=1e-7)
                for key in ('edges', 'observed', 'observed_cov', 'chi2_mask'):
                    saved[root+key] = pt[0][root+key]
                for name, source, key, label, _, _ in METHODS:
                    dataset = pt if source == 'pt' else mass
                    for suffix in ('', '_cov', '_closure_bins', '_residual_cov'):
                        saved[root+name+suffix] = dataset[0][root+key+suffix]
                    row = dict(dataset[1][(cat, tag, ch, 'VR', cal, 'mass', key)])
                    obs, obsvar = saved[root+'observed']
                    pred, predvar = saved[root+name]
                    residual_cov = saved[root+name+'_residual_cov']
                    chi2, rank = recompute_statistics(obs, pred, residual_cov,
                                                     saved[root+'chi2_mask'].astype(bool))
                    np.testing.assert_allclose(chi2, float(row['chi2']), rtol=1e-10, atol=1e-8)
                    assert rank == int(row['rank'])
                    np.testing.assert_allclose(np.diag(saved[root+name+'_cov']), predvar, rtol=1e-10, atol=1e-8)
                    crossdiag = (obsvar+predvar-np.diag(residual_cov))/2
                    good = pred > 0
                    joint_ratio_var = (obsvar[good]/pred[good]**2+
                        obs[good]**2*predvar[good]/pred[good]**4-
                        2*obs[good]*crossdiag[good]/pred[good]**3)
                    np.testing.assert_allclose(saved[root+name+'_closure_bins'][0, good],
                                               obs[good]/pred[good], rtol=1e-10, atol=1e-8)
                    np.testing.assert_allclose(saved[root+name+'_closure_bins'][1, good],
                                               joint_ratio_var, rtol=1e-9, atol=1e-8)
                    row.update(chi2=chi2, rank=rank,
                        chi2_pvalue=float(chi2_distribution.sf(chi2, rank)) if rank else np.nan)
                    row.update(method=name, method_label=label, source_study=source)
                    rows.append(row)
    return saved, rows


def plot(data, summaries, cat, tag, ch, cal, compact=False, axes=None):
    if axes is None:
        fig, axes = plt.subplots(2, 1, figsize=(11, 10), sharex=True,
                                gridspec_kw={'height_ratios': [3, 1.15], 'hspace': .05})
        fig.subplots_adjust(top=.91, bottom=.15, left=.14, right=.97)
    else:
        fig = axes[0].figure
    ax, lower = axes
    root = prefix(cat, tag, ch, cal)
    x = data[root+'edges'].copy()
    if np.isinf(x[-1]):
        x[-1] = x[-2]+max(250, min(800, (x[-2]-x[0])*.35))
    centers, widths = (x[1:]+x[:-1])/2, np.diff(x)
    observation, variance = data[root+'observed']
    error = np.sqrt(np.maximum(0., variance))
    ax.errorbar(centers, observation, yerr=error, fmt='o', color='black',
                ms=3 if compact else 4, label='Data - TTplusWW - VZ', zorder=8)
    ymin, ymax = min(0., np.min(observation-error)), np.max(observation+error)
    visible = []
    step = lambda v: np.r_[v, v[-1]]
    for i, (name, _, _, label, color, style) in enumerate(METHODS):
        value, var = data[root+name]
        sigma = np.sqrt(np.maximum(0., var))
        row = summaries[(cat, tag, ch, cal, name)]
        full_label = label if compact else label+f"  (Z/pred. = {float(row['closure']):.3f})"
        ax.step(x, step(value), where='post', color=color, lw=1.9,
                ls=style, label=full_label, zorder=3+i*.2)
        ax.fill_between(x, step(value-sigma), step(value+sigma), step='post',
                        color=color, alpha=.075, lw=0)
        ratio, ratio_var = data[root+name+'_closure_bins']
        valid = data[root+'chi2_mask'].astype(bool) & np.isfinite(ratio) & np.isfinite(ratio_var)
        lower.errorbar(centers[valid]+(i-2)*.045*widths[valid], ratio[valid],
                       yerr=np.sqrt(np.maximum(0., ratio_var[valid])), fmt='o',
                       color=color, ms=2 if compact else 3, elinewidth=.7, capsize=0)
        visible.extend(ratio[valid].tolist())
        ymin, ymax = min(ymin, np.min(value-sigma)), max(ymax, np.max(value+sigma))
    ax.set_ylim(ymin-.035*(ymax-ymin), max(1., ymax)*(1.8 if compact else 2.05))
    hep.cms.label('Work in progress', data=True, lumi=59.8, com=13,
                  fontsize=10 if compact else 19, ax=ax)
    ax.text(.03, .97, f'2018, VR\n{cat.capitalize()}, {LABELS[tag]}\n{CH[ch]}',
            transform=ax.transAxes, va='top', fontsize=10 if compact else 14)
    if not compact:
        ax.legend(loc='upper right', bbox_to_anchor=(.99, .98), fontsize=11, frameon=False)
    lower.axhline(1., color='.35', lw=.8)
    lower.set_ylim(min(.7, min(visible)-.1) if visible else .7,
                   max(1.3, max(visible)+.1) if visible else 1.3)
    lower.set_ylabel('Obs. / pred.', fontsize=11 if compact else 17)
    candidates = (300, 500, 750, 1100, 1600, 2300)
    ticks = [v for v in candidates if x[0] <= v < x[-2]]+[centers[-1]]
    lower.set_xticks(ticks)
    lower.set_xticklabels([f'{v:g}' for v in ticks[:-1]]+[rf'$\geq {x[-2]:g}$'])
    lower.set_xlabel(r'$m_{2\ell2q}$ [GeV]' if cat == 'resolved' else r'$m_{2\ell J}$ [GeV]',
                     fontsize=13 if compact else 20)
    ax.set_ylabel('Events / bin', fontsize=12 if compact else 19)
    ax.tick_params(labelbottom=False)
    for item in axes:
        item.set_xlim(x[0], x[-1])
        item.tick_params(which='both', direction='in', top=True, right=True,
                         labelsize=9 if compact else 13)
        item.grid(False)
    if not compact:
        qualifier = ('CR + VR correction calibration: VR partly in-sample.' if cal == 'inclVR'
                     else 'CR-only correction calibration: VR held out of fitted factor values.')
        fig.text(.5, .075, qualifier+' No ZH subtraction.', ha='center', fontsize=10)
        fig.text(.5, .045, 'Bands: prediction stat. errors; ratio errors include observation/prediction covariance.', ha='center', fontsize=10)
        fig.text(.5, .018, 'Alpha source: CR excluding VR. Final bin includes overflow; original mass bins retained.', ha='center', fontsize=10)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pt', type=Path, required=True)
    parser.add_argument('--mass', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a fresh comparison output folder')
    pt, mass = load(args.pt), load(args.mass)
    data, rows = prepare(pt, mass)
    summaries = {(r['category'], r['tag'], r['channel'], r['calibration'], r['method']): r for r in rows}
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output/'results.npz', **data)
    with (args.output/'summary.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plt.style.use(hep.style.CMS)
    plt.rcParams.update({'font.size': 14, 'axes.labelpad': 7, 'pdf.fonttype': 42, 'axes.linewidth': 1.1})
    def save(fig, stem, bundle):
        fig.savefig(args.output/(stem+'.pdf'), bbox_inches='tight')
        fig.savefig(args.output/(stem+'.png'), dpi=180, bbox_inches='tight')
        bundle.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        print('Saved '+stem, flush=True)
    for cal in CALIBRATIONS:
      with PdfPages(args.output/f'all_VR_method_comparison_{cal}.pdf') as bundle:
        for ch in CHANNELS:
            fig = plt.figure(figsize=(18, 12))
            outer = fig.add_gridspec(2, 3, wspace=.32, hspace=.38)
            for i, (cat, tag) in enumerate((c, t) for c in ('resolved', 'merged') for t in TAGS):
                grid = outer[i//3, i%3].subgridspec(2, 1, height_ratios=[3, 1.15], hspace=.05)
                ax = fig.add_subplot(grid[0])
                lower = fig.add_subplot(grid[1], sharex=ax)
                plot(data, summaries, cat, tag, ch, cal, compact=True, axes=(ax, lower))
            handles, labels = fig.axes[0].get_legend_handles_labels()
            fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .995),
                       ncol=3, fontsize=11, frameon=False)
            qualifier = ('CR+VR correction calibration: VR partly in-sample' if cal == 'inclVR'
                         else 'CR-only correction calibration: VR held out of fitted factor values')
            fig.text(.5, .018, qualifier+'; alpha source excludes VR; stat. only; no ZH subtraction.', ha='center', fontsize=11)
            fig.subplots_adjust(top=.875, bottom=.10, left=.07, right=.985)
            save(fig, f'{ch}_VR_methods_summary_{cal}', bundle)
        for cat in ('resolved', 'merged'):
            for tag in TAGS:
                for ch in CHANNELS:
                    save(plot(data, summaries, cat, tag, ch, cal), f'{cat}_{tag}_{ch}_VR_methods_{cal}', bundle)
    for ch in CHANNELS:
      for cal in CALIBRATIONS:
        table = [f'# VR chi2 comparison, {ch}, {cal}', '',
            'Entries are chi2/rank, with full statistical residual covariance and fixed original mass bins.', '',
            '| Category | Direct DY | DY + pT | DY + mass | Nominal alpha | Alpha + pT |',
            '|---|---:|---:|---:|---:|---:|']
        for cat in ('resolved', 'merged'):
          for tag in TAGS:
            values = []
            for name, *_ in METHODS:
                row = summaries[(cat, tag, ch, cal, name)]
                values.append(f"{float(row['chi2']):.2f}/{int(row['rank'])}")
            table.append('| '+cat+' '+LABELS[tag]+' | '+' | '.join(values)+' |')
        (args.output/f'{ch}_VR_chi2_{cal}.md').write_text('\n'.join(table)+'\n')
    metadata = dict(date='2026-10-07', source_sha256=pt[2]['source_sha256'],
        calibrations={'inclVR': 'C+VR factors: VR partly in-sample',
                      'heldoutVR': 'C-only fitted factor values; VR held out'},
        alpha_source='C only, disjoint from VR', ZH_subtracted=False,
        original_mass_grid_preserved=True, independent_channel_category_corrections=True,
        uncertainty='Copied full prediction covariance and joint observation/prediction ratio errors; methods are correlated, not independent measurements',
        input_sha256={name: hashlib.sha256((folder/'results.npz').read_bytes()).hexdigest()
                      for name, folder in [('pt', args.pt), ('mass', args.mass)]},
        input_folders={'pt': str(args.pt.resolve()), 'mass': str(args.mass.resolve())})
    for name, folder in [('pt', args.pt), ('mass', args.mass)]:
        for file in ('metadata.json', 'summary.csv'):
            shutil.copy2(folder/file, args.output/(name+'_'+file))
    shutil.copy2(Path(__file__), args.output/Path(__file__).name)
    shutil.copy2(Path(__file__).with_name('VR_METHOD_COMPARISON.md'), args.output/'README.md')
    (args.output/'metadata.json').write_text(json.dumps(metadata, indent=2))
    assert len(list(args.output.glob('*.png'))) == 42
    assert len(list(args.output.glob('*.pdf'))) == 44
    assert len(rows) == 180
    validation = dict(status='passed', matched_VR_category_channels_calibrations=36,
                      checked_fields_per_category_channel=8, recomputed_chi2=180,
                      validated_joint_ratio_errors=180, predictions=180,
                      individual_png_pdf_pairs=42, bundle_pages={'inclVR': 21, 'heldoutVR': 21},
                      no_refit=True, source_sha256=pt[2]['source_sha256'])
    (args.output/'validation.json').write_text(json.dumps(validation, indent=2))
    with zipfile.ZipFile(args.output.with_suffix('.zip'), 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(args.output.iterdir()):
            if path.is_file():
                archive.write(path, arcname=args.output.name+'/'+path.name)
    print(json.dumps(validation), flush=True)


if __name__ == '__main__':
    main()
