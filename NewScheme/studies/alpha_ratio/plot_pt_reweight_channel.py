"""CMS-style comparisons for independent category/channel DY pT corrections."""
import argparse
import csv
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
from pt_reweight_core import TAGS
from pt_reweight_channel import prefix, CHANNELS
from plot_pt_reweight import display_edges, step, LABELS, CH, COLORS


def comparison(data, summaries, cat, tag, ch, target, cal, axis, method, compact=False, axes=None):
    if axes is None:
        fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True,
                                gridspec_kw={'height_ratios': [3, 1], 'hspace': .05})
        fig.subplots_adjust(top=.91, bottom=.16, left=.15, right=.96)
    else:
        fig = axes[0].figure
    ax, lower = axes
    root = prefix(cat, tag, ch, target, cal, axis)
    edges = data[root+'edges']
    x = display_edges(edges)
    center = (x[:-1]+x[1:])/2
    observation, obsvar = data[root+'observed']
    ax.errorbar(center, observation, yerr=np.sqrt(obsvar), fmt='o', color='black',
                ms=3 if compact else 4, label='Data - TTplusWW - VZ', zorder=5)
    nominal = 'dy_nominal' if method == 'direct' else 'alpha_nominal'
    corrected = 'dy_reweighted' if method == 'direct' else 'alpha_reweighted'
    curves = [('nominal', nominal, 'Nominal DY' if method == 'direct' else 'Nominal alpha'),
              ('category', corrected, 'pT reweighted DY' if method == 'direct' else 'Reweighted alpha')]
    if method == 'direct' and target != 'calibration':
        curves.append(('normalization', 'dy_normonly', 'Normalization only'))
    ymax = np.max(observation+np.sqrt(obsvar))
    ymin = min(0., np.min(observation-np.sqrt(obsvar)))
    visible = []
    for color_key, key, label in curves:
        values, variance = data[root+key]
        error = np.sqrt(np.maximum(0., variance))
        color = COLORS[color_key]
        ax.step(x, step(values), where='post', color=color, lw=1.8,
                ls='--' if color_key != 'nominal' else '-', label=label)
        if color_key != 'normalization':
            ax.fill_between(x, step(values-error), step(values+error), step='post',
                            color=color, alpha=.13, lw=0)
            ratio, ratio_variance = data[root+key+'_closure_bins']
            mask = data[root+'chi2_mask'].astype(bool) & np.isfinite(ratio)
            offset = -.06 if color_key == 'nominal' else .06
            lower.errorbar(center[mask]+offset*np.diff(x)[mask], ratio[mask],
                           yerr=np.sqrt(np.maximum(0., ratio_variance[mask])), fmt='o',
                           ms=2.5 if compact else 3, color=color, elinewidth=.8)
            visible.extend(ratio[mask].tolist())
        ymax = max(ymax, np.max(values+error))
        ymin = min(ymin, np.min(values-error))
    ax.set_ylim(ymin-.035*(ymax-ymin), max(1., ymax)*1.6)
    hep.cms.label('Work in progress', data=True, lumi=59.8, com=13,
                  fontsize=10 if compact else 19, ax=ax)
    region = ('CR + VR' if cal == 'inclVR' else 'CR only') if target == 'calibration' else target
    ax.text(.035, .965, f'2018, {region}\n{cat.capitalize()}, {LABELS[tag]}\n{CH[ch]}',
            transform=ax.transAxes, va='top', fontsize=10 if compact else 14)
    nominal_row = summaries[(cat, tag, ch, target, cal, axis, nominal)]
    corrected_row = summaries[(cat, tag, ch, target, cal, axis, corrected)]
    ax.text(.975, .965, 'Observed/predicted:\n'
            f"Nominal: {float(nominal_row['closure']):.3f}\n"
            f"Reweighted: {float(corrected_row['closure']):.3f}",
            transform=ax.transAxes, ha='right', va='top', fontsize=8 if compact else 11)
    if not compact:
        ax.legend(loc='upper right', bbox_to_anchor=(.99, .76), fontsize=11, frameon=False)
    lower.axhline(1., color='.3', lw=.9)
    lower.set_ylabel('Obs. / pred.', fontsize=11 if compact else 17)
    lower.set_ylim(min(.7, min(visible)-.1) if visible else .7,
                   max(1.3, max(visible)+.1) if visible else 1.3)
    candidates = (300, 500, 750, 1100, 1600, 2300) if axis == 'mass' else (0, 200, 300, 400, 500, 650, 850)
    ticks = [v for v in candidates if x[0] <= v < x[-2]]+[center[-1]]
    lower.set_xticks(ticks)
    lower.set_xticklabels([f'{v:g}' for v in ticks[:-1]]+[rf'$\geq {x[-2]:g}$'])
    xlabel = (r'$m_{2\ell2q}$ [GeV]' if cat == 'resolved' else r'$m_{2\ell J}$ [GeV]') if axis == 'mass' else (
              r'$p_T^{qq}$ [GeV]' if cat == 'resolved' else r'$p_T^J$ [GeV]')
    lower.set_xlabel(xlabel, fontsize=13 if compact else 20)
    ax.set_ylabel('Events / bin', fontsize=12 if compact else 19)
    ax.tick_params(labelbottom=False)
    for item in axes:
        item.set_xlim(x[0], x[-1])
        item.tick_params(which='both', direction='in', top=True, right=True, labelsize=9 if compact else 13)
        item.grid(False)
    if not compact:
        qualifier = ('In-sample calibration: exact closure by construction.' if target == 'calibration'
                     else 'VR enters pT calibration: in-sample diagnostic' if target == 'VR' and cal == 'inclVR'
                     else 'SR and VR held out of calibration' if cal == 'heldoutVR' else 'SR held out of calibration')
        fig.text(.5, .085, 'Independent pT correction per category and lepton channel', ha='center', fontsize=11)
        fig.text(.5, .055, qualifier, ha='center', fontsize=10)
        fig.text(.5, .02, 'Statistical errors include fit/data covariance; last bin includes overflow.', ha='center', fontsize=10)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    args = parser.parse_args()
    output = args.input
    if (output/'all_pt_reweight_plots.pdf').exists():
        raise FileExistsError('Plot bundle already exists')
    with np.load(output/'results.npz') as archive:
        data = {key: archive[key] for key in archive.files}
    with (output/'summary.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    summaries = {tuple(row[k] for k in ('category', 'tag', 'channel', 'target', 'calibration', 'axis', 'prediction')): row for row in rows}
    plt.style.use(hep.style.CMS)
    plt.rcParams.update({'font.size': 14, 'axes.labelpad': 7, 'pdf.fonttype': 42, 'axes.linewidth': 1.1})
    counts = {}

    def save(fig, stem, bundle, kind):
        fig.savefig(output/(stem+'.pdf'), bbox_inches='tight')
        fig.savefig(output/(stem+'.png'), dpi=180, bbox_inches='tight')
        bundle.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        counts[kind] = counts.get(kind, 0)+1

    with PdfPages(output/'all_pt_reweight_plots.pdf') as bundle:
        for ch in CHANNELS:
            for method, target in [('direct', 'SR'), ('alpha', 'SR'), ('alpha', 'VR')]:
                fig = plt.figure(figsize=(17, 11.5))
                outer = fig.add_gridspec(2, 3, wspace=.32, hspace=.38)
                for i, (cat, tag) in enumerate((c, t) for c in ('resolved', 'merged') for t in TAGS):
                    grid = outer[i//3, i%3].subgridspec(2, 1, height_ratios=[3, 1], hspace=.05)
                    ax = fig.add_subplot(grid[0])
                    lower = fig.add_subplot(grid[1], sharex=ax)
                    comparison(data, summaries, cat, tag, ch, target, 'inclVR', 'mass', method,
                               compact=True, axes=(ax, lower))
                handles, labels = fig.axes[0].get_legend_handles_labels()
                fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .995),
                           ncol=4 if method == 'direct' else 3, fontsize=11, frameon=False)
                qualifier = 'VR enters calibration' if target == 'VR' else 'SR held out of calibration'
                fig.text(.5, .015, f'Channel-dependent, category-specific correction; {qualifier}. Statistical errors only.',
                         ha='center', fontsize=11)
                fig.subplots_adjust(top=.90, bottom=.10, left=.075, right=.985)
                save(fig, f'{ch}_{target}_{method}_summary', bundle, 'mass_and_factors')
        for cat in ('resolved', 'merged'):
            for tag in TAGS:
                for ch in CHANNELS:
                    for method, target in [('direct', 'SR'), ('alpha', 'SR'), ('alpha', 'VR')]:
                        fig = comparison(data, summaries, cat, tag, ch, target, 'inclVR', 'mass', method)
                        save(fig, f'{cat}_{tag}_{ch}_{target}_{method}', bundle, 'mass_and_factors')
                    fig = comparison(data, summaries, cat, tag, ch, 'VR', 'heldoutVR', 'mass', 'alpha')
                    save(fig, f'{cat}_{tag}_{ch}_VR_alpha_heldoutVR', bundle, 'mass_and_factors')
                for ch in ('ee', 'mumu'):
                    fig, ax = plt.subplots(figsize=(10, 7.8))
                    fig.subplots_adjust(left=.15, right=.96, top=.91, bottom=.16)
                    ranges = []
                    for cal, color, label in [('inclVR', '#2166ac', 'CR + VR calibration'),
                                              ('heldoutVR', '#b2182b', 'CR-only calibration')]:
                        root = f'{cat}__{tag}__{ch}__SF__{cal}__'
                        x = display_edges(data[root+'edges'])
                        values, variance = data[root+'values']
                        error = np.sqrt(variance)
                        ax.step(x, step(values), where='post', color=color, lw=1.8, label=label)
                        ax.fill_between(x, step(values-error), step(values+error), step='post', color=color, alpha=.16, lw=0)
                        ranges.append(x[-1])
                    hep.cms.label('Work in progress', data=True, lumi=59.8, com=13, fontsize=19, ax=ax)
                    ax.text(.035, .97, f'2018, {cat.capitalize()}\n{LABELS[tag]}\n{CH[ch]}',
                            transform=ax.transAxes, va='top', fontsize=14)
                    ax.legend(loc='upper right', bbox_to_anchor=(.99, .80), fontsize=12)
                    ax.set_xlabel(r'$p_T^{qq}$ [GeV]' if cat == 'resolved' else r'$p_T^J$ [GeV]', fontsize=20)
                    ax.set_ylabel('(Data - TTplusWW - VZ) / DY', fontsize=17)
                    ax.axhline(1., color='.3', lw=.8)
                    ax.set_xlim(0., max(ranges))
                    ymin, ymax = ax.get_ylim()
                    ax.set_ylim(min(0., ymin), ymax*1.5)
                    ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=13)
                    fig.text(.5, .035, 'Independent channel/category fit; statistical errors only; final bin includes overflow.',
                             ha='center', fontsize=10)
                    save(fig, f'{cat}_{tag}_{ch}_pt_scale_factors', bundle, 'mass_and_factors')
    with PdfPages(output/'all_SR_pt_reweight.pdf') as bundle:
        for cat in ('resolved', 'merged'):
            for tag in TAGS:
                for ch in CHANNELS:
                    fig = comparison(data, summaries, cat, tag, ch, 'SR', 'inclVR', 'pt', 'direct')
                    save(fig, f'{cat}_{tag}_{ch}_SR_direct_pt', bundle, 'SR_pt')
    with PdfPages(output/'all_CR_pt_reweight.pdf') as bundle:
        for cal in ('inclVR', 'heldoutVR'):
            for cat in ('resolved', 'merged'):
                for tag in TAGS:
                    for ch in CHANNELS:
                        fig = comparison(data, summaries, cat, tag, ch, 'calibration', cal, 'pt', 'direct')
                        save(fig, f'{cat}_{tag}_{ch}_{cal}_direct_pt', bundle, 'CR_pt')
    assert counts == {'mass_and_factors': 93, 'SR_pt': 18, 'CR_pt': 36}, counts
    for name in ('pt_reweight_core.py', 'plot_pt_reweight.py', 'pt_reweight_channel.py',
                 'plot_pt_reweight_channel.py', 'test_pt_reweight.py', 'test_pt_reweight_channel.py'):
        shutil.copy2(Path(__file__).parent/name, output/name)
    document = Path(__file__).parent/'PT_CHANNEL_STUDY.md'
    if not document.exists():
        document = Path(__file__).parent/'README.md'
    shutil.copy2(document, output/'README.md')
    metadata = json.loads((output/'metadata.json').read_text())
    metadata['individual_png_pdf_pairs'] = sum(counts.values())
    metadata['bundle_pages'] = counts
    metadata['script_sha256'] = {name: hashlib.sha256((output/name).read_bytes()).hexdigest()
                                 for name in metadata.get('script_sha256', [])}
    for name in ('pt_reweight_channel.py', 'plot_pt_reweight_channel.py', 'pt_reweight_core.py'):
        metadata['script_sha256'][name] = hashlib.sha256((output/name).read_bytes()).hexdigest()
    (output/'metadata.json').write_text(json.dumps(metadata, indent=2))
    pngs = list(output.glob('*.png'))
    assert len(pngs) == 147 and all(p.with_suffix('.pdf').exists() for p in pngs)
    with zipfile.ZipFile(output.with_suffix('.zip'), 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(output.iterdir()):
            if path.is_file():
                archive.write(path, arcname=output.name+'/'+path.name)
    print('Saved 147 individual PNG/PDF pairs, 93/18/36-page bundles, and ZIP.', flush=True)


if __name__ == '__main__':
    main()
