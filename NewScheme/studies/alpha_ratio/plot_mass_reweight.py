"""CMS plots and self-contained package for the 2018 CR mass-reweighting test."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
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
from pt_reweight_channel import CHANNELS, CALIBRATIONS
from plot_pt_reweight import display_edges, step, LABELS, CH
from plot_pt_reweight_channel import comparison


def massplot(data, summaries, cat, tag, ch, target, method, cal, compact=False, axes=None):
    fig = comparison(data, summaries, cat, tag, ch, target, cal, 'mass', method, compact, axes)
    for ax in fig.axes:
        for line in ax.lines:
            label = line.get_label()
            if label == 'pT reweighted DY':
                line.set_label('Mass-reweighted DY')
            elif label == 'Reweighted alpha':
                line.set_label('Mass-reweighted alpha')
        legend = ax.get_legend()
        if legend is not None:
            for text in legend.get_texts():
                text.set_text(text.get_text().replace('pT reweighted DY', 'Mass-reweighted DY')
                              .replace('Reweighted alpha', 'Mass-reweighted alpha'))
    if not compact:
        for text in fig.texts:
            text.set_text(text.get_text().replace('Independent pT correction', 'Independent mass correction')
                          .replace('VR enters pT calibration', 'VR enters mass calibration')
                          .replace('In-sample calibration: exact closure by construction.', 'In-sample calibration: exact closure only in the coarser correction bins.'))
        if method == 'alpha':
            fig.text(.5, .12, 'Nominal and mass-reweighted alpha coincide (exact binwise cancellation).', ha='center', fontsize=10)
        elif target == 'SR' and cal == 'inclVR':
            fig.text(.5, .12, 'Coarse CR mass correction applied on the preserved comparison grid.', ha='center', fontsize=10)
    return fig


def finalize(folder):
    assert len(list(folder.glob('*.png'))) == 135
    assert len(list(folder.glob('*.pdf'))) == 137
    assert all(p.with_suffix('.pdf').exists() for p in folder.glob('*.png'))
    metadata = json.loads((folder/'metadata.json').read_text())
    assert hashlib.sha256((folder/'source_histograms.npz').read_bytes()).hexdigest() == metadata['source_sha256']
    with np.load(folder/'results.npz') as data:
        for key in data.files:
            if key.endswith('__alpha_reweighted'):
                nominal = key[:-len('alpha_reweighted')]+'alpha_nominal'
                np.testing.assert_allclose(data[key], data[nominal], rtol=1e-10, atol=1e-7)
                np.testing.assert_allclose(data[key+'_cov'], data[nominal+'_cov'], rtol=1e-10, atol=1e-7)
    result = subprocess.run([sys.executable, '-m', 'unittest', 'test_pt_reweight', 'test_mass_reweight', '-v'],
                            cwd=str(folder.resolve()), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, check=True)
    (folder/'tests.log').write_text(result.stdout)
    doc = Path(__file__).parent/'MASS_REWEIGHT_STUDY.md'
    if doc.exists():
        shutil.copy2(doc, folder/'README.md')
    if Path(__file__).parent.resolve() != folder.resolve():
        shutil.copy2(Path(__file__), folder/Path(__file__).name)
    checks = metadata['checks']
    validation = dict(status='passed', regression_tests=15, individual_png_pdf_pairs=135,
        PDF_files=137, bundle_pages={'comparisons_and_factors': 93, 'calibration': 42},
        calibration_closures=sum(c['test'] == 'mass_calibration_closure' for c in checks),
        alpha_cancellation_checks=sum(c['test'] == 'alpha_factor_cancellation' for c in checks),
        direct_equals_alpha_checks=sum(c['test'] == 'direct_equals_nominal_alpha' for c in checks),
        source_sha256=metadata['source_sha256'])
    (folder/'validation.json').write_text(json.dumps(validation, indent=2))
    metadata['script_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.glob('*.py')}
    (folder/'metadata.json').write_text(json.dumps(metadata, indent=2))
    with zipfile.ZipFile(folder.with_suffix('.zip'), 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(folder.iterdir()):
            if path.is_file():
                archive.write(path, arcname=folder.name+'/'+path.name)
    print(json.dumps(validation), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--finalize-only', action='store_true')
    args = parser.parse_args()
    output = args.input
    if args.finalize_only:
        finalize(output)
        return
    if (output/'all_mass_reweight_plots.pdf').exists():
        raise FileExistsError('Plots already exist')
    with np.load(output/'results.npz') as archive:
        data = {key: archive[key] for key in archive.files}
    with (output/'summary.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    summaries = {tuple(row[k] for k in ('category', 'tag', 'channel', 'target', 'calibration', 'axis', 'prediction')): row for row in rows}
    plt.style.use(hep.style.CMS)
    plt.rcParams.update({'font.size': 14, 'axes.labelpad': 7, 'pdf.fonttype': 42, 'axes.linewidth': 1.1})
    count = 0
    def save(fig, stem, bundle):
        nonlocal count
        fig.savefig(output/(stem+'.pdf'), bbox_inches='tight')
        fig.savefig(output/(stem+'.png'), dpi=180, bbox_inches='tight')
        bundle.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        count += 1

    def overview(ch, target, method, cal):
        fig = plt.figure(figsize=(17, 11.5))
        outer = fig.add_gridspec(2, 3, wspace=.32, hspace=.38)
        for i, (cat, tag) in enumerate((c, t) for c in ('resolved', 'merged') for t in TAGS):
            grid = outer[i//3, i%3].subgridspec(2, 1, height_ratios=[3, 1], hspace=.05)
            ax = fig.add_subplot(grid[0])
            lower = fig.add_subplot(grid[1], sharex=ax)
            massplot(data, summaries, cat, tag, ch, target, method, cal, compact=True, axes=(ax, lower))
        handles, labels = fig.axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .995),
                   ncol=4 if method == 'direct' else 3, fontsize=11, frameon=False)
        if method == 'alpha':
            note = 'Nominal and mass-reweighted alpha coincide exactly; target VR enters primary calibration.' if target == 'VR' else 'Nominal and mass-reweighted alpha coincide exactly; SR held out of calibration.'
        elif target == 'calibration':
            note = 'In-sample mass calibration: closure imposed only in the coarser correction bins.'
        else:
            note = 'Coarse CR mass correction on preserved comparison grid; SR held out of calibration.'
        fig.text(.5, .015, note+' Statistical errors only.', ha='center', fontsize=11)
        fig.subplots_adjust(top=.90, bottom=.10, left=.075, right=.985)
        return fig

    with PdfPages(output/'all_mass_reweight_plots.pdf') as bundle:
        for ch in CHANNELS:
            for method, target in [('direct', 'SR'), ('alpha', 'SR'), ('alpha', 'VR')]:
                save(overview(ch, target, method, 'inclVR'), f'{ch}_{target}_{method}_summary', bundle)
        for cat in ('resolved', 'merged'):
            for tag in TAGS:
                for ch in CHANNELS:
                    for method, target in [('direct', 'SR'), ('alpha', 'SR'), ('alpha', 'VR')]:
                        save(massplot(data, summaries, cat, tag, ch, target, method, 'inclVR'), f'{cat}_{tag}_{ch}_{target}_{method}', bundle)
                    save(massplot(data, summaries, cat, tag, ch, 'VR', 'alpha', 'heldoutVR'), f'{cat}_{tag}_{ch}_VR_alpha_heldoutVR', bundle)
                for ch in ('ee', 'mumu'):
                    fig, ax = plt.subplots(figsize=(10, 7.8))
                    fig.subplots_adjust(left=.15, right=.96, top=.91, bottom=.16)
                    x = None
                    for cal, color, label in [('inclVR', '#2166ac', 'CR + VR calibration'), ('heldoutVR', '#b2182b', 'CR-only calibration')]:
                        root = f'{cat}__{tag}__{ch}__SF__{cal}__'
                        x = display_edges(data[root+'edges'])
                        val, var = data[root+'values']
                        error = np.sqrt(var)
                        ax.step(x, step(val), where='post', color=color, lw=1.8, label=label)
                        ax.fill_between(x, step(val-error), step(val+error), step='post', color=color, alpha=.16, lw=0)
                    hep.cms.label('Work in progress', data=True, lumi=59.8, com=13, fontsize=19, ax=ax)
                    ax.text(.035, .97, f'2018, {cat.capitalize()}\n{LABELS[tag]}\n{CH[ch]}', transform=ax.transAxes, va='top', fontsize=14)
                    ax.legend(loc='upper right', bbox_to_anchor=(.99, .80), fontsize=12)
                    ax.set_xlabel(r'$m_{2\ell2q}$ [GeV]' if cat == 'resolved' else r'$m_{2\ell J}$ [GeV]', fontsize=20)
                    ax.set_ylabel('(Data - TTplusWW - VZ) / DY', fontsize=17)
                    ax.axhline(1., color='.3', lw=.8)
                    ax.set_xlim(x[0], x[-1])
                    ymin, ymax = ax.get_ylim()
                    ax.set_ylim(min(0., ymin), ymax*1.5)
                    center = (x[:-1]+x[1:])/2
                    ticks = [v for v in (300, 500, 750, 1100, 1600) if v < x[-2]]+[center[-1]]
                    ax.set_xticks(ticks)
                    ax.set_xticklabels([f'{v:g}' for v in ticks[:-1]]+[rf'$\geq {x[-2]:g}$'])
                    ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=13)
                    fig.text(.5, .035, 'Independent category/channel mass correction; stat. errors only; final bin includes overflow.', ha='center', fontsize=10)
                    save(fig, f'{cat}_{tag}_{ch}_mass_scale_factors', bundle)
    assert count == 93
    with PdfPages(output/'all_CR_mass_reweight.pdf') as bundle:
        for cal in CALIBRATIONS:
            for ch in CHANNELS:
                save(overview(ch, 'calibration', 'direct', cal), f'{ch}_{cal}_CR_mass_summary', bundle)
            for cat in ('resolved', 'merged'):
                for tag in TAGS:
                    for ch in CHANNELS:
                        save(massplot(data, summaries, cat, tag, ch, 'calibration', 'direct', cal), f'{cat}_{tag}_{ch}_{cal}_CR_mass', bundle)
    assert count == 135
    for name in ('mass_reweight_dy.py', 'plot_mass_reweight.py', 'test_mass_reweight.py',
                 'pt_reweight_core.py', 'pt_reweight_channel.py', 'plot_pt_reweight_channel.py',
                 'plot_pt_reweight.py', 'test_pt_reweight.py'):
        shutil.copy2(Path(__file__).parent/name, output/name)
    doc = Path(__file__).parent/'MASS_REWEIGHT_STUDY.md'
    if not doc.exists():
        doc = Path(__file__).parent/'README.md'
    shutil.copy2(doc, output/'README.md')
    finalize(output)


if __name__ == '__main__':
    main()
