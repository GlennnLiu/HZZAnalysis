"""Compare five saved SR predictions, without refitting or changing mass bins."""
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

TAGS = ('btag', 'untag', 'vbftag')
CHANNELS = ('ee', 'mumu', 'combined')
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


def prefix(cat, tag, ch):
    return '__'.join((cat, tag, ch, 'SR', 'inclVR', 'mass'))+'__'


def prepare(pt, mass):
    """Fail closed if nominal inputs, bins or covariance differ across studies."""
    if pt[2]['source_sha256'] != mass[2]['source_sha256']:
        raise ValueError('Different primitive histogram inputs')
    saved, rows = {}, []
    for cat in ('resolved', 'merged'):
        for tag in TAGS:
            for ch in CHANNELS:
                root = prefix(cat, tag, ch)
                for key in ('edges', 'observed', 'observed_cov', 'chi2_mask',
                            'dy_nominal', 'dy_nominal_cov', 'alpha_nominal', 'alpha_nominal_cov'):
                    np.testing.assert_allclose(pt[0][root+key], mass[0][root+key], rtol=1e-10, atol=1e-7)
                for key in ('edges', 'observed', 'observed_cov', 'chi2_mask'):
                    saved[root+key] = pt[0][root+key]
                for name, source, key, label, _, _ in METHODS:
                    dataset = pt if source == 'pt' else mass
                    for suffix in ('', '_cov', '_closure_bins'):
                        saved[root+name+suffix] = dataset[0][root+key+suffix]
                    row = dict(dataset[1][(cat, tag, ch, 'SR', 'inclVR', 'mass', key)])
                    row.update(method=name, method_label=label, source_study=source)
                    rows.append(row)
    return saved, rows


def plot(data, summaries, cat, tag, ch, compact=False, axes=None):
    if axes is None:
        fig, axes = plt.subplots(2, 1, figsize=(11, 10), sharex=True,
                                gridspec_kw={'height_ratios': [3, 1.15], 'hspace': .05})
        fig.subplots_adjust(top=.91, bottom=.15, left=.14, right=.97)
    else:
        fig = axes[0].figure
    ax, lower = axes
    root = prefix(cat, tag, ch)
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
        row = summaries[(cat, tag, ch, name)]
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
    ax.text(.03, .97, f'2018, SR\n{cat.capitalize()}, {LABELS[tag]}\n{CH[ch]}',
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
        fig.text(.5, .075, 'CR + VR calibration; independent category/channel corrections; no ZH subtraction.', ha='center', fontsize=10)
        fig.text(.5, .045, 'Bands: prediction stat. errors; ratio errors include observation/prediction covariance.', ha='center', fontsize=10)
        fig.text(.5, .018, 'SR held out; final bin includes overflow. Near-overlapping curves are retained.', ha='center', fontsize=10)
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
    summaries = {(r['category'], r['tag'], r['channel'], r['method']): r for r in rows}
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
    with PdfPages(args.output/'all_SR_method_comparison.pdf') as bundle:
        for ch in CHANNELS:
            fig = plt.figure(figsize=(18, 12))
            outer = fig.add_gridspec(2, 3, wspace=.32, hspace=.38)
            for i, (cat, tag) in enumerate((c, t) for c in ('resolved', 'merged') for t in TAGS):
                grid = outer[i//3, i%3].subgridspec(2, 1, height_ratios=[3, 1.15], hspace=.05)
                ax = fig.add_subplot(grid[0])
                lower = fig.add_subplot(grid[1], sharex=ax)
                plot(data, summaries, cat, tag, ch, compact=True, axes=(ax, lower))
            handles, labels = fig.axes[0].get_legend_handles_labels()
            fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .995),
                       ncol=3, fontsize=11, frameon=False)
            fig.text(.5, .018, 'CR+VR calibration; channel/category dependent; SR held out; statistical errors only; no ZH subtraction.', ha='center', fontsize=11)
            fig.subplots_adjust(top=.875, bottom=.10, left=.07, right=.985)
            save(fig, ch+'_SR_methods_summary', bundle)
        for cat in ('resolved', 'merged'):
            for tag in TAGS:
                for ch in CHANNELS:
                    save(plot(data, summaries, cat, tag, ch), f'{cat}_{tag}_{ch}_SR_methods', bundle)
    metadata = dict(date='2026-10-07', source_sha256=pt[2]['source_sha256'],
        calibration='C+VR for reweighting and alpha SR source', ZH_subtracted=False,
        original_mass_grid_preserved=True, independent_channel_category_corrections=True,
        uncertainty='Copied full prediction covariance and joint observation/prediction ratio errors; methods are correlated, not independent measurements',
        input_sha256={name: hashlib.sha256((folder/'results.npz').read_bytes()).hexdigest()
                      for name, folder in [('pt', args.pt), ('mass', args.mass)]},
        input_folders={'pt': str(args.pt.resolve()), 'mass': str(args.mass.resolve())})
    for name, folder in [('pt', args.pt), ('mass', args.mass)]:
        for file in ('metadata.json', 'summary.csv'):
            shutil.copy2(folder/file, args.output/(name+'_'+file))
    shutil.copy2(Path(__file__), args.output/Path(__file__).name)
    shutil.copy2(Path(__file__).with_name('SR_METHOD_COMPARISON.md'), args.output/'README.md')
    (args.output/'metadata.json').write_text(json.dumps(metadata, indent=2))
    assert len(list(args.output.glob('*.png'))) == 21
    assert len(list(args.output.glob('*.pdf'))) == 22
    assert len(rows) == 90
    validation = dict(status='passed', matched_SR_category_channels=18,
                      checked_fields_per_category_channel=8, predictions=90,
                      individual_png_pdf_pairs=21, bundle_pages=21,
                      no_refit=True, source_sha256=pt[2]['source_sha256'])
    (args.output/'validation.json').write_text(json.dumps(validation, indent=2))
    with zipfile.ZipFile(args.output.with_suffix('.zip'), 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(args.output.iterdir()):
            if path.is_file():
                archive.write(path, arcname=args.output.name+'/'+path.name)
    print(json.dumps(validation), flush=True)


if __name__ == '__main__':
    main()
