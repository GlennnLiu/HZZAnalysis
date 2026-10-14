#!/usr/bin/env python3
"""Plot saved SR DY predictions only; never read SR observation arrays."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import mplhep as hep


CATEGORIES = [(cat, tag) for cat in ('resolved', 'merged')
              for tag in ('btag', 'untag', 'vbftag')]
CHANNELS = ('ee', 'mumu', 'combined')
LABELS = {'btag': 'b-tag', 'untag': 'untagged', 'vbftag': 'VBF-tag'}
CHANNEL_LABELS = {'ee': r'$ee$', 'mumu': r'$\mu\mu$',
                  'combined': r'$ee+\mu\mu$'}
BLUE, RED, BAND = '#2166ac', '#b2182b', '#a6c8e3'


def load_predictions(archive):
    records = {}
    # Only these four arrays per category are read, not any SR observation.
    with np.load(archive, allow_pickle=False) as saved:
        for channel in CHANNELS:
            for cat, tag in CATEGORIES:
                prefix = '__'.join((cat, tag, channel, 'SR')) + '__'
                edges = saved[prefix + 'edges'].copy()
                old = saved[prefix + 'pred_nozh'].copy()
                new = saved[prefix + 'pred_withzh'].copy()
                transferred = saved[prefix + 'transferred_zh'].copy()
                assert old.shape == new.shape == transferred.shape == (2, len(edges)-1)
                assert edges[0] == 300 and np.isposinf(edges[-1])
                assert np.all(np.diff(edges) > 0)
                assert np.all(np.isfinite(old)) and np.all(np.isfinite(new))
                assert np.all(old[1] >= 0) and np.all(new[1] >= 0)
                np.testing.assert_allclose(new[0], old[0]-transferred[0], rtol=1e-10, atol=1e-10)
                sigma = np.sqrt(old[1])
                delta = new[0]-old[0]
                assert np.all(delta[sigma == 0] == 0), 'Nonzero shift with zero uncertainty'
                scaled = np.divide(delta, sigma, out=np.zeros_like(delta), where=sigma > 0)
                display = edges.copy()
                display[-1] = edges[-2] + max(350., min(1000., (edges[-2]-300.)*.35))
                records[cat, tag, channel] = dict(edges=edges, display=display,
                                                old=old, new=new, sigma=sigma, scaled=scaled)
    return records


def step(values):
    return np.r_[values, values[-1]]


def draw_pair(ax, lower, record, cat, tag, channel, compact=False):
    x, old, new, sigma = record['display'], record['old'][0], record['new'][0], record['sigma']
    ax.fill_between(x, step(old-sigma), step(old+sigma), step='post',
                    color=BAND, alpha=.65, linewidth=0, label=r'Stat. uncertainty ($\pm1\sigma$)', zorder=1)
    ax.step(x, step(old), where='post', color=BLUE, linewidth=1.8,
            label='Without ZH subtraction', zorder=2)
    ax.step(x, step(new), where='post', color=RED, linewidth=1.7, linestyle='--',
            label='With ZH subtraction in CR', zorder=3)
    ax.axhline(0, color='0.5', linewidth=.6, zorder=0)
    ax.set_ylabel('Events / bin', fontsize=12 if compact else 19)
    hep.cms.label('Work in progress', data=True, lumi=59.8, com=13,
                  loc=0, fontsize=11 if compact else 19, ax=ax)
    ax.text(.035, .965, f'2018, SR\n{cat.capitalize()}, {LABELS[tag]}\n{CHANNEL_LABELS[channel]}',
            transform=ax.transAxes, ha='left', va='top', fontsize=10 if compact else 14)
    ymax = max(float(np.max(old+sigma)), float(np.max(new)), 1.)
    ymin = min(float(np.min(old-sigma)), float(np.min(new)), 0.)
    ax.set_ylim(ymin-.035*(ymax-ymin), ymax*1.62)
    old_total, new_total = float(old.sum()), float(new.sum())
    stat_total = float(np.sqrt(record['old'][1].sum()))
    annotation = f'Yield (without / with):\n{old_total:,.2f} ± {stat_total:,.2f} (stat)\n{new_total:,.2f}'
    ax.text(.98, .97, annotation, transform=ax.transAxes, ha='right', va='top',
            fontsize=8 if compact else 11)

    lower.axhspan(-1, 1, color=BAND, alpha=.65, linewidth=0)
    lower.axhline(0, color=BLUE, linewidth=.9)
    lower.step(x, step(record['scaled']), where='post', color=RED, linestyle='--', linewidth=1.7)
    centers = (x[:-1]+x[1:])/2
    lower.plot(centers, record['scaled'], 'o', color=RED, markersize=3)
    lower.set_ylim(-1.12, 1.12)
    lower.set_yticks([-1, 0, 1])
    lower.set_ylabel(r'$\Delta/\sigma_{\mathrm{stat}}$', fontsize=12 if compact else 19)
    lower.text(.98, .95, r'$\max\,|\Delta|/\sigma_{\mathrm{stat}} = $' + f"{np.max(np.abs(record['scaled'])):.3f}",
               transform=lower.transAxes, ha='right', va='top', fontsize=9 if compact else 12)
    masslabel = r'$m_{2\ell2q}$ [GeV]' if cat == 'resolved' else r'$m_{2\ell J}$ [GeV]'
    lower.set_xlabel(masslabel, fontsize=13 if compact else 20)
    # Sparse ticks avoid crowded low-mass labels in the six-category overview.
    ticks = [v for v in (300, 500, 750, 1100, 1600, 2300) if v < x[-2]] + [centers[-1]]
    labels = [f'{v:g}' for v in ticks[:-1]] + [rf'$\geq {x[-2]:g}$']
    lower.set_xticks(ticks)
    lower.set_xticklabels(labels, fontsize=9 if compact else 13)
    lower.text(.015, .06, 'Last bin includes overflow', transform=lower.transAxes,
               fontsize=8 if compact else 10, ha='left', va='bottom')
    for axis in (ax, lower):
        axis.set_xlim(x[0], x[-1])
        axis.grid(False)
        axis.tick_params(which='both', direction='in', top=True, right=True,
                         labelsize=9 if compact else 13)
        axis.tick_params(which='major', length=5 if compact else 7)
        axis.tick_params(which='minor', length=2.5 if compact else 3.5)
    ax.tick_params(labelbottom=False)


def save(fig, output, stem, bundle):
    fig.savefig(output / (stem + '.pdf'), bbox_inches='tight')
    fig.savefig(output / (stem + '.png'), dpi=180, bbox_inches='tight')
    bundle.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError(f'Refusing to overwrite nonempty output: {args.output}')
    records = load_predictions(args.input)
    args.output.mkdir(parents=True, exist_ok=True)
    plt.style.use(hep.style.CMS)
    plt.rcParams.update({'font.size': 14, 'axes.unicode_minus': False,
                         'axes.labelpad': 7, 'axes.linewidth': 1.1,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    rows = []
    with PdfPages(args.output / 'all_SR_ZH_statband.pdf') as bundle:
        # First page: all six categories with combined ee+mu mu predictions.
        fig = plt.figure(figsize=(17, 11.5))
        outer = fig.add_gridspec(2, 3, wspace=.32, hspace=.38)
        for i, (cat, tag) in enumerate(CATEGORIES):
            grid = outer[i//3, i%3].subgridspec(2, 1, height_ratios=[3, 1], hspace=.04)
            ax = fig.add_subplot(grid[0])
            lower = fig.add_subplot(grid[1], sharex=ax)
            draw_pair(ax, lower, records[cat, tag, 'combined'], cat, tag, 'combined', compact=True)
        handles, labels = fig.axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .995), ncol=3, frameon=False, fontsize=13)
        fig.text(.5, .015, r'DY alpha-method estimate; statistical uncertainties only.  $\Delta =$ with $-$ without.  No SR data plotted.',
                 ha='center', fontsize=12)
        fig.subplots_adjust(top=.91, bottom=.10, left=.075, right=.985)
        save(fig, args.output, 'combined_SR_ZH_statband_summary', bundle)

        for channel in CHANNELS:
            for cat, tag in CATEGORIES:
                record = records[cat, tag, channel]
                fig, (ax, lower) = plt.subplots(2, 1, figsize=(10, 9), sharex=True,
                    gridspec_kw=dict(height_ratios=[3, 1], hspace=.05))
                draw_pair(ax, lower, record, cat, tag, channel)
                ax.legend(loc='upper right', bbox_to_anchor=(.99, .80), fontsize=12, frameon=False)
                fig.text(.5, .025, r'DY alpha-method estimate; statistical uncertainties only.  $\Delta =$ with $-$ without.', ha='center', fontsize=11)
                fig.subplots_adjust(top=.92, bottom=.15, left=.15, right=.96)
                stem = f'{cat}_{tag}_{channel}_SR_ZH_statband'
                save(fig, args.output, stem, bundle)
                old_yield = float(record['old'][0].sum())
                new_yield = float(record['new'][0].sum())
                sigma_yield = float(np.sqrt(record['old'][1].sum()))
                rows.append(dict(category=cat, tag=tag, channel=channel,
                    yield_without=old_yield, yield_without_stat=sigma_yield,
                    yield_with=new_yield, yield_with_stat=float(np.sqrt(record['new'][1].sum())),
                    delta=new_yield-old_yield,
                    absolute_yield_shift_over_stat=abs(new_yield-old_yield)/sigma_yield,
                    max_absolute_bin_shift_over_stat=float(np.max(np.abs(record['scaled']))),
                    bins=len(record['scaled']), bins_outside_band=int(np.count_nonzero(np.abs(record['scaled'])>1)),
                    file=stem+'.pdf'))
    with (args.output/'summary.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    provenance = dict(source=str(args.input.resolve()),
        source_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        read_arrays=['edges', 'pred_nozh', 'pred_withzh', 'transferred_zh'],
        observation_arrays_read=False, target='SR', channels=list(CHANNELS),
        alpha='Nominal DY SR/CR, as saved in the original broad-CR study',
        plots=19, pdf_bundle_pages=19, style='mplhep.style.CMS',
        mplhep_version=hep.__version__, cms_label='Work in progress',
        high_mass_region_shading=False, luminosity_label_fb=59.8,
        maximum_bin_shift_sigma=max(row['max_absolute_bin_shift_over_stat'] for row in rows))
    (args.output/'metadata.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (args.output/'README.md').write_text('''# 2018 SR: ZH subtraction compared with statistical uncertainty

Start with `combined_SR_ZH_statband_summary.pdf` (six categories, ee + mumu).
Individual PDFs and PNGs cover resolved/merged × btag/untag/vbftag × ee/mumu/combined.
`all_SR_ZH_statband.pdf` collects all 19 figures in a single PDF.

Blue: P0 = alpha * (Data_CR - TTplusWW_CR - VZ_CR).
Red dashed: P1 = alpha * (Data_CR - TTplusWW_CR - VZ_CR - ZHbb_CR).
Light-blue band: P0 +/- sqrt(Var(P0)), the statistical uncertainty of the original estimate.
The band includes finite CR data, subtracted MC and DY alpha statistics; no systematics.
Lower panels show (P1 - P0) / sqrt(Var(P0)); the shaded interval is [-1, +1].
This is a comparison to the ORIGINAL prediction precision, not a significance test of the difference.
Before/after predictions share the same data and alpha and are strongly correlated.

Only DY estimates are shown: ZH in the target SR is NOT added to the red curve.
No SR Data-TT-VZ points or observation arrays are read by this plotting script.
The subtraction is the available ZH_HToBB_ZToLL sample, not every Higgs background.
The broad hadronic-mass CR includes [105,135]; this is not an m2l2q control region.
All search masses >=300 GeV are retained. The final bin labelled >= its lower edge
includes the entire open-ended high-mass tail; its finite plotted width is only a
display convention. No high-mass shading or hatching is applied: no region is excluded
or blinded. Shading is used only for the statistical uncertainty band.
Signed bin contents are preserved. ee+mumu predictions and variances are sums of
the separately transferred channels, not a combined alpha applied to combined CR data.

`summary.csv` records yields, yield statistical errors and maximum bin shifts in sigma units.
`metadata.json` records the source archive and its SHA-256 digest.
`plot_sr_zh_statband.py` reproduces these plots from the saved histogram archive:

    python plot_sr_zh_statband.py --input ../zh_subtraction_2018_broadCR/histograms.npz --output NEW_OUTPUT_DIRECTORY

The output directory must be empty. No reduced trees or production alpha files are modified.

Presentation: mplhep.style.CMS, also used in HZZAnalysis/FitTemplate/make_1Dtemplate.py.
CMS Work in progress marks an internal diagnostic, not an approved public result.
The label is NOT Simulation: these are estimates derived from CR data, with MC transfer
and subtraction. Luminosity is shown as 59.8 fb^-1 (rounded from the unchanged 59.83).
Official mplhep styling documentation: https://mplhep.readthedocs.io/en/latest/guide_styling/
''')
    shutil.copy2(Path(__file__), args.output/'plot_sr_zh_statband.py')
    archive_path = args.output.with_suffix('.zip')
    if archive_path.exists():
        raise FileExistsError(f'Refusing to overwrite archive: {archive_path}')
    with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(args.output.iterdir()):
            if path.is_file():
                archive.write(path, arcname=str(Path(args.output.name)/path.name))
    print(json.dumps(dict(output=str(args.output), zip=str(archive_path), **provenance), indent=2), flush=True)


if __name__ == '__main__':
    main()
