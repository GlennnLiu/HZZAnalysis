"""Plot calibration-region DY pT closure from the existing source cache."""
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
from pt_reweight_core import Space, Estimate, TAGS, fit_pt, compare
from plot_pt_reweight import COLORS, LABELS, CH, display_edges, step


def aggregate(item, groups):
    return Estimate(np.array([item.value[g].sum() for g in groups]),
                    np.stack([item.jac[g].sum(axis=0) for g in groups]))


def weighted(fine_dy, fit, groups):
    lookup = np.empty(len(fine_dy.value), dtype=int)
    for i, group in enumerate(fit['groups']):
        lookup[group] = i
    expanded = Estimate(fit['factor'].value[lookup], fit['factor'].jac[lookup])
    return aggregate(fine_dy * expanded, groups)


def plot(cat, tag, ch, cal, edges, obs, curves, variance, summaries):
    fig, (ax, lower) = plt.subplots(2, 1, figsize=(10, 9), sharex=True,
        gridspec_kw={'height_ratios': [3, 1], 'hspace': .05})
    fig.subplots_adjust(top=.91, bottom=.18, left=.15, right=.96)
    x = display_edges(edges)
    center = (x[:-1] + x[1:]) / 2
    error = np.sqrt(np.diag(obs.covariance(variance)))
    ax.errorbar(center, obs.value, yerr=error, fmt='o', color='black', ms=4,
                label='Data - TTplusWW - VZ', zorder=5)
    ymax = np.max(obs.value + error)
    ymin = min(0., np.min(obs.value - error))
    visible = []
    for index, (key, prediction, label) in enumerate(curves):
        color = COLORS[key]
        err = np.sqrt(np.maximum(0., np.diag(prediction.covariance(variance))))
        ax.step(x, step(prediction.value), where='post', color=color, lw=1.8,
                ls='--' if key == 'category' else '-', label=label)
        ax.fill_between(x, step(prediction.value-err), step(prediction.value+err),
                        step='post', color=color, alpha=.13, lw=0)
        ymax = max(ymax, np.max(prediction.value + err))
        ymin = min(ymin, np.min(prediction.value - err))
        good = prediction.value > 0
        ratio = Estimate(obs.value[good], obs.jac[good]) / Estimate(
            prediction.value[good], prediction.jac[good])
        ratio_err = np.sqrt(np.maximum(0., np.diag(ratio.covariance(variance))))
        offset = (index-(len(curves)-1)/2) * .10
        lower.errorbar(center[good] + offset*np.diff(x)[good], ratio.value,
                       yerr=ratio_err, fmt='o', ms=3, color=color, elinewidth=.8)
        visible.extend(ratio.value.tolist())
    ax.set_ylim(ymin-.035*(ymax-ymin), max(1., ymax)*1.6)
    hep.cms.label('Work in progress', data=True, lumi=59.8, com=13, fontsize=19, ax=ax)
    title_tag = 'All three tags pooled' if tag == 'shared' else LABELS[tag]
    region = 'CR + VR calibration' if cal == 'inclVR' else 'CR-only calibration'
    ax.text(.035, .965, f'2018, {region}\n{cat.capitalize()}, {title_tag}\n{CH[ch]}',
            transform=ax.transAxes, va='top', fontsize=14)
    main_key = 'shared' if tag == 'shared' else 'category'
    ax.text(.975, .965, 'Observed/predicted:\n'
            f"Nominal: {summaries['nominal']['closure']:.3f}\n"
            f"After correction: {summaries[main_key]['closure']:.3f}",
            transform=ax.transAxes, ha='right', va='top', fontsize=11)
    ax.legend(loc='upper right', bbox_to_anchor=(.99, .75), fontsize=11, frameon=False)
    lower.axhline(1., color='.3', lw=.9)
    lower.set_ylim(min(.7, min(visible)-.1), max(1.3, max(visible)+.1))
    lower.set_ylabel('Obs. / pred.', fontsize=17)
    ticks = [v for v in (0, 200, 300, 400, 500, 650, 850) if v < x[-2]] + [center[-1]]
    lower.set_xticks(ticks)
    lower.set_xticklabels([f'{v:g}' for v in ticks[:-1]] + [rf'$\geq {x[-2]:g}$'])
    lower.set_xlabel(r'$p_T^{qq}$ [GeV]' if cat == 'resolved' else r'$p_T^J$ [GeV]', fontsize=20)
    ax.set_ylabel('Events / bin', fontsize=19)
    for axis in (ax, lower):
        axis.set_xlim(x[0], x[-1])
        axis.tick_params(which='both', direction='in', top=True, right=True, labelsize=13)
        axis.grid(False)
    ax.tick_params(labelbottom=False)
    fig.text(.5, .095, 'In-sample calibration: bin-by-bin closure enforced for the fitted sample.',
             ha='center', fontsize=10)
    fig.text(.5, .06, 'Factors fitted to ee + '+r'$\mu\mu$'+
             '; separate channels need not close exactly.', ha='center', fontsize=10)
    fig.text(.5, .025, 'Statistical errors include shared fit/data covariance; last bin includes overflow.',
             ha='center', fontsize=10)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    source_file = args.source / 'source_histograms.npz'
    with np.load(source_file) as archive:
        arrays = {key: archive[key] for key in archive.files}
    plt.style.use(hep.style.CMS)
    plt.rcParams.update({'font.size': 14, 'pdf.fonttype': 42, 'axes.linewidth': 1.1})
    rows, bins, saved, assertions = [], [], {}, []
    with PdfPages(args.output / 'all_CR_pt_reweight.pdf') as bundle:
        for cal, regions in [('inclVR', [0, 2]), ('CRonly', [0])]:
            for cat in ('resolved', 'merged'):
                space = Space(arrays, cat)
                fits = {tag: fit_pt(space, list(range(3)) if tag == 'shared' else [TAGS.index(tag)], regions)
                        for tag in ('shared',)+TAGS}
                finegroups = [np.array([i]) for i in range(space.shape[-1])]
                for tag in ('shared',)+TAGS:
                    tags = list(range(3)) if tag == 'shared' else [TAGS.index(tag)]
                    fit = fits[tag]
                    groups, edges = fit['groups'], fit['edges']
                    for ch, channels in [('combined', [0, 1]), ('ee', [0]), ('mumu', [1])]:
                        obs = space.linear(0, tags, channels, regions, ptgroups=groups)
                        fine_dy = space.linear(1, tags, channels, regions, ptgroups=finegroups)
                        nominal = aggregate(fine_dy, groups)
                        own = weighted(fine_dy, fit, groups)
                        main_key = 'shared' if tag == 'shared' else 'category'
                        curves = [('nominal', nominal, 'Nominal DY')]
                        if tag != 'shared':
                            curves.append(('shared', weighted(fine_dy, fits['shared'], groups),
                                           'pT reweight: shared tags'))
                        curves.append((main_key, own, 'pT reweight: pooled tags' if tag == 'shared'
                                       else 'pT reweight: per category'))
                        if ch == 'combined':
                            np.testing.assert_allclose(own.value, obs.value, rtol=1e-12, atol=1e-9)
                            np.testing.assert_allclose(own.jac, obs.jac, rtol=0, atol=1e-12)
                            assertions.append(dict(category=cat, tag=tag, calibration=cal,
                                max_absolute_yield_residual=float(np.max(np.abs(own.value-obs.value))),
                                max_absolute_jacobian_residual=float(np.max(np.abs(own.jac-obs.jac)))))
                        root = '__'.join((cat, tag, ch, cal))+'__'
                        saved[root+'edges'] = edges
                        saved[root+'observed'] = np.stack([obs.value, np.diag(obs.covariance(space.variance))])
                        summaries = {}
                        for key, item, label in curves:
                            result = compare(obs, item, space.variance)
                            summaries[key] = result
                            rows.append(dict(category=cat, tag=tag, channel=ch, calibration=cal,
                                             prediction=key, **result))
                            saved[root+key] = np.stack([item.value, np.diag(item.covariance(space.variance))])
                            saved[root+key+'_cov'] = item.covariance(space.variance)
                            residual = obs-item
                            saved[root+key+'_residual_cov'] = residual.covariance(space.variance)
                            good = item.value > 0
                            rr = Estimate(obs.value[good], obs.jac[good])/Estimate(item.value[good], item.jac[good])
                            ratio = np.full((2, len(groups)), np.nan)
                            ratio[0, good] = rr.value
                            ratio[1, good] = np.diag(rr.covariance(space.variance))
                            saved[root+key+'_closure_bins'] = ratio
                            for i in range(len(groups)):
                                bins.append(dict(category=cat, tag=tag, channel=ch, calibration=cal,
                                    prediction=key, pt_low=edges[i], pt_high=edges[i+1],
                                    observed=obs.value[i], predicted=item.value[i],
                                    closure=ratio[0, i], closure_stat=np.sqrt(max(0., ratio[1, i]))))
                        fig = plot(cat, tag, ch, cal, edges, obs, curves, space.variance, summaries)
                        stem = f'{cat}_{tag}_{ch}_{cal}_direct_pt'
                        fig.savefig(args.output/(stem+'.pdf'), bbox_inches='tight')
                        fig.savefig(args.output/(stem+'.png'), dpi=180, bbox_inches='tight')
                        bundle.savefig(fig, bbox_inches='tight')
                        plt.close(fig)
    for name, content in [('cr_pt_summary.csv', rows), ('cr_pt_bins.csv', bins)]:
        with (args.output/name).open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(content[0]))
            writer.writeheader()
            writer.writerows(content)
    np.savez_compressed(args.output/'cr_pt_distributions.npz', **saved)
    metadata = dict(source=str(source_file.resolve()), source_sha256=hashlib.sha256(source_file.read_bytes()).hexdigest(),
                    individual_png_pdf_pairs=48, bundle_pages=48, exact_closure_checks=assertions,
                    fit_channels='ee+mumu', calibration_regions={'inclVR': 'C+VR', 'CRonly': 'C'},
                    zero_fitted_sample_ratio_error='Shared calibration/data covariance cancels; not independent validation.')
    (args.output/'validation.json').write_text(json.dumps(metadata, indent=2))
    for filename in ('plot_cr_pt_reweight.py', 'pt_reweight_core.py', 'plot_pt_reweight.py'):
        shutil.copy2(Path(__file__).parent/filename, args.output/filename)
    shutil.copy2(Path(__file__).parent/'CR_PT_PLOTS.md', args.output/'README.md')
    with zipfile.ZipFile(args.output.with_suffix('.zip'), 'w', zipfile.ZIP_DEFLATED) as archive:
        for filename in sorted(args.output.iterdir()):
            archive.write(filename, arcname=args.output.name+'/'+filename.name)
    print('Saved 48 PNG/PDF pairs and 48-page bundle. All 16 combined fitted samples close in values and Jacobians.', flush=True)


if __name__ == '__main__':
    main()
