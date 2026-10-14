"""2018 broad-sideband DY alpha closure, before/after ZHbb subtraction.

Uses freshly recomputed nominal coarse-bin DY ratios, not stored production
alpha histograms. All search masses >=300, including the high-mass tail, enter.
Combined-channel predictions are sums of the independently transferred ee/mumu
predictions. Original inputs and shared analysis configuration remain unchanged.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
for name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[name] = '1'
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CHANNELS = ['ee', 'mumu']
CATEGORIES = ['resolved', 'merged']
TAGS = ['btag', 'untag', 'vbftag']
PROCESSES = ['Data', 'DY', 'TTplusWW', 'VZ', 'ZHbb']
INITIAL_EDGES = {
    'resolved': [300, 400, 500, 600, 750, 900, 1100, 1400, 1800, 2300, 3500, np.inf],
    'merged': [300, 500, 750, 1100, 1600, 3500, np.inf]}


def hist(a, edges, weight='weight', extra=None):
    mask = np.ones(len(a['mass']), dtype=bool) if extra is None else extra
    m, w = a['mass'][mask], a[weight][mask]
    return np.stack([np.histogram(m, edges, weights=w)[0],
                     np.histogram(m, edges, weights=w*w)[0]])


def divide(n, d):
    return np.divide(n, d, out=np.full(np.broadcast(n, d).shape, np.nan), where=d > 0)


def derive(h):
    dc, vdc = h['CR', 'DY']; dt, vdt = h['target', 'DY']
    if np.any(dc <= 0) or np.any(dt <= 0):
        raise ValueError('Use positive DY-only merged bins; do not invent alpha in empty bins')
    c = h['CR', 'Data'][0]-h['CR', 'TTplusWW'][0]-h['CR', 'VZ'][0]
    vc = h['CR', 'Data'][1]+h['CR', 'TTplusWW'][1]+h['CR', 'VZ'][1]
    t = h['target', 'Data'][0]-h['target', 'TTplusWW'][0]-h['target', 'VZ'][0]
    vt = h['target', 'Data'][1]+h['target', 'TTplusWW'][1]+h['target', 'VZ'][1]
    hc, vhc = h['CR', 'ZHbb']; ht, vht = h['target', 'ZHbb']
    hu, vhu = h['CR', 'ZHbb_unit_tag']; htu, vhtu = h['target', 'ZHbb_unit_tag']
    alpha = dt/dc
    va = vdt/(dc*dc)+dt*dt*vdc/(dc**4)
    out = dict(h)
    def add(name, value, variance): out[name] = np.stack([value, variance])
    add('alpha', alpha, va)
    add('cr_nozh', c, vc); add('cr_withzh', c-hc, vc+vhc)
    add('obs_nozh', t, vt); add('obs_withzh', t-ht, vt+vht)
    add('pred_nozh', alpha*c, alpha*alpha*vc+c*c*va)
    add('pred_withzh', alpha*(c-hc), alpha*alpha*(vc+vhc)+(c-hc)**2*va)
    add('pred_total_withzh', alpha*(c-hc)+ht, out['pred_withzh'][1]+vht)
    add('transferred_zh', alpha*hc, alpha*alpha*vhc+hc*hc*va)
    add('transferred_zh_unit_tag', alpha*hu, alpha*alpha*vhu+hu*hu*va)
    add('pred_withzh_unit_tag', alpha*(c-hu), alpha*alpha*(vc+vhu)+(c-hu)**2*va)
    add('obs_withzh_unit_tag', t-htu, vt+vhtu)
    add('delta_dy_prediction', -alpha*hc, out['transferred_zh'][1])
    add('delta_total_prediction', ht-alpha*hc, vht+out['transferred_zh'][1])
    add('residual_nozh', t-alpha*c, vt+out['pred_nozh'][1])
    add('residual_withzh', t-ht-alpha*(c-hc), vt+vht+out['pred_withzh'][1])
    # The differences share all baseline data/other-background fluctuations;
    # they cancel. Do not add the before/after uncertainties as independent.
    add('delta_closure_residual', alpha*hc-ht, vht+out['transferred_zh'][1])
    out['cov_before_after_prediction'] = alpha*alpha*vc+c*(c-hc)*va
    out['cov_transferred_zh_pred_nozh'] = hc*c*va
    return out


def combine(results):
    result = {key: sum(r[key] for r in results) for key in results[0]}
    dc, vdc = result['CR', 'DY']; dt, vdt = result['target', 'DY']
    # Diagnostic combined DY ratio only; predictions above are sums of the
    # channel-specific transfers, not this ratio times the combined CR.
    result['alpha'] = np.stack([dt/dc, vdt/dc**2+dt**2*vdc/dc**4])
    return result


def read_inputs(args, cat, tag, audit):
    arrays = {}
    for ch in CHANNELS:
        for reg in ['CR', 'SR', 'VR']:
            for proc in PROCESSES:
                base = args.zh_base if proc == 'ZHbb' else (args.cr_base if reg == 'CR' else args.base)
                path = base/f'{proc}_{ch}_{cat}_{tag}_2018_{reg}.root'
                fields = ['mass', 'mass2q', 'weight']+(['weight_unit_tag'] if proc == 'ZHbb' else [])
                with uproot.open(path) as f: a = f['tree'].arrays(fields, library='np')
                if not all(np.all(np.isfinite(v)) for v in a.values()):
                    raise ValueError('Nonfinite reduced values: %s' % path)
                q = a['mass2q']; m = a['mass']
                sr = (q > 70) & (q < 105); vr = (q > 135) & (q < 165)
                if reg == 'SR' and not np.all(sr): raise ValueError('Wrong SR reduction: %s' % path)
                if reg == 'VR' and not np.all(vr): raise ValueError('Wrong VR reduction: %s' % path)
                if reg == 'CR' and np.any(sr): raise ValueError('CR includes SR: %s' % path)
                selected = (m >= 300) & (q > 40) & (q < 180)
                if reg == 'CR': selected &= ((q < 70) | (q > 105)) & ~vr
                audit.append(dict(path=str(path), bytes=path.stat().st_size,
                    mtime=path.stat().st_mtime, entries=len(q), selected=int(selected.sum()),
                    selected_overflow=int((selected & (m >= 3500)).sum()),
                    selected_h_window=int((selected & (q > 105) & (q <= 135)).sum()),
                    vr_removed=int((vr & (m >= 300)).sum()) if reg == 'CR' else 0,
                    selected_weight=float(a['weight'][selected].sum())))
                arrays[ch, reg, proc] = {k: v[selected].astype(float) for k, v in a.items()}
    return arrays


def adaptive_edges(arrays, cat):
    edges = list(INITIAL_EDGES[cat])
    merges = []
    while True:
        bad = np.zeros(len(edges)-1, dtype=bool)
        for ch in CHANNELS:
            for reg in ['CR', 'SR', 'VR']:
                y, _ = hist(arrays[ch, reg, 'DY'], edges)
                bad |= y <= 0
        if not bad.any(): return np.array(edges), merges
        if len(edges) == 2: raise ValueError('Nonpositive integrated DY yield')
        index = int(np.flatnonzero(bad)[-1])
        remove = index if index == len(edges)-2 else index+1
        merges.append(float(edges[remove]))
        edges.pop(remove)


def make_histograms(arrays, ch, target, edges):
    out = {}
    for reg, key in [('CR', 'CR'), (target, 'target')]:
        for proc in PROCESSES:
            a = arrays[ch, reg, proc]
            out[key, proc] = hist(a, edges)
            if proc == 'ZHbb':
                out[key, 'ZHbb_unit_tag'] = hist(a, edges, weight='weight_unit_tag')
                if reg == 'CR':
                    out[key, 'ZHbb_Hwindow'] = hist(a, edges, extra=(a['mass2q'] > 105) & (a['mass2q'] <= 135))
    return out


def yield_row(r, cat, tag, ch, target):
    row = dict(category=cat, tag=tag, channel=ch, target=target)
    for key, value in r.items():
        if isinstance(key, tuple): name = '_'.join(key)
        else: name = key
        if isinstance(key, str) and key.startswith('cov_'): continue
        if key == 'alpha': continue
        row[name] = float(value[0].sum())
        row[name+'_stat'] = float(np.sqrt(value[1].sum()))
    row['fraction_dy_prediction_removed'] = float(divide(row['transferred_zh'], row['pred_nozh']))
    row['fraction_target_zh'] = float(divide(row['target_ZHbb'], row['obs_nozh']))
    row['fraction_total_prediction_change'] = float(divide(row['delta_total_prediction'], row['pred_nozh']))
    row['closure_nozh'] = float(divide(row['obs_nozh'], row['pred_nozh']))
    row['closure_withzh'] = float(divide(row['obs_withzh'], row['pred_withzh']))
    row['closure_total_withzh'] = float(divide(row['obs_nozh'], row['pred_total_withzh']))
    for suffix in ['nozh', 'withzh']:
        v = r['residual_'+suffix][1]
        row['chi2_'+suffix] = float(np.sum(r['residual_'+suffix][0]**2/v))
    row['diagnostic_bins'] = len(r['obs_nozh'][0])
    return row


def save_figure(fig, output, stem):
    for ext in ['pdf', 'png']: fig.savefig(output/(stem+'.'+ext), dpi=170)
    plt.close(fig)


def distribution_plot(r, edges, cat, tag, ch, target, output):
    displayed = edges.copy()
    displayed[-1] = displayed[-2] + max(350., min(1000., (displayed[-2]-300)*.35))
    x = .5*(displayed[:-1]+displayed[1:]); dx = np.diff(displayed)
    fig, (ax, ar, ad) = plt.subplots(3, 1, figsize=(10, 10), sharex=True,
                    gridspec_kw={'height_ratios': [3.4, 1.5, 1.2]}, constrained_layout=True)
    no, vn = r['obs_nozh']; yes, vy = r['obs_withzh']
    ax.errorbar(x-.025*dx, no, np.sqrt(vn), fmt='o', color='black', ms=4,
                label='Target Data − TT/WW/tZq − VZ')
    ax.errorbar(x+.025*dx, yes, np.sqrt(vy), fmt='o', mfc='white', color='tab:green', ms=4,
                label='Target Data − TT/WW/tZq − VZ − ZHbb')
    for key, color, label, ls in [
        ('pred_nozh', 'tab:blue', 'Alpha estimate: no ZH subtraction', '-'),
        ('pred_withzh', 'tab:green', 'Alpha estimate: ZH subtracted in CR', '-'),
        ('pred_total_withzh', 'tab:red', 'Corrected DY estimate + target ZHbb', '--'),
        (('target', 'DY'), 'tab:orange', 'Direct DY MC (reference)', ':')]:
        y, v = r[key]
        ax.stairs(y, displayed, color=color, label=label, ls=ls, lw=1.7)
        if key in ['pred_nozh', 'pred_withzh']:
            ax.fill_between(displayed, np.r_[y-np.sqrt(v), (y-np.sqrt(v))[-1]],
                np.r_[y+np.sqrt(v), (y+np.sqrt(v))[-1]], step='post', color=color, alpha=.12)
    ax.set_title(f'2018 {cat}, {tag}, {ch} — {target}\nBroad hadronic sideband retains 105–135 GeV; full search mass ≥300 GeV')
    ax.set_ylabel('Events / bin')
    ax.legend(fontsize=8, ncol=2, loc='upper right')
    ax.axhline(0, color='gray', lw=.6)
    lo = min(0., float(np.min(no-np.sqrt(vn))), float(np.min(yes-np.sqrt(vy))))
    ax.set_ylim(bottom=lo*1.15 if lo < 0 else 0)
    for obs, pred, color, offset, label in [
        ('obs_nozh', 'pred_nozh', 'tab:blue', -.03, 'Before: target / alpha'),
        ('obs_withzh', 'pred_withzh', 'tab:green', .03, 'After: DY-only target / alpha')]:
        y, v = r[obs]; p, vp = r[pred]
        ar.errorbar(x+offset*dx, divide(y, p), divide(np.sqrt(v), p), fmt='o', color=color, ms=4, label=label)
        band = divide(np.sqrt(vp), p)
        ar.fill_between(displayed, np.r_[1-band, (1-band)[-1]],
                        np.r_[1+band, (1+band)[-1]], step='post', color=color, alpha=.12)
    ar.axhline(1, color='gray', ls='--'); ar.set_ylabel('Target / estimate')
    ar.set_ylim(-.5, 2.5); ar.legend(fontsize=8, loc='upper right')
    p, vp = r['pred_nozh']; xfer, vxfer = r['transferred_zh']; ht, vht = r['target', 'ZHbb']
    cov = r['cov_transferred_zh_pred_nozh']
    for value, variance, covariance, color, label in [
        (-xfer, vxfer, -cov, 'tab:green', 'DY prediction change'),
        (ht-xfer, vht+vxfer, -cov, 'tab:red', 'DY+ZH prediction change')]:
        frac = divide(value, p)
        vf = divide(variance, p*p)+divide(value*value*vp, p**4)-divide(2*value*covariance, p**3)
        ad.errorbar(x, 100*frac, 100*np.sqrt(np.maximum(vf, 0)), fmt='o-', ms=3, color=color, label=label)
    ad.axhline(0, color='gray', ls='--'); ad.set_ylabel('Change [%]')
    ad.legend(fontsize=8, loc='upper right')
    ad.set_xlabel('m(2ℓ2q) [GeV]; last bin includes the entire high-mass tail')
    for axis in [ax, ar, ad]:
        axis.axvspan(displayed[-2], displayed[-1], color='gray', alpha=.06)
    ticks = list(displayed[:-1]); labels = [f'{v:g}' for v in ticks]
    if len(ticks) > 9: ticks, labels = ticks[::2]+[displayed[-2]], labels[::2]+[f'{displayed[-2]:g}']
    ad.set_xticks(ticks, labels, rotation=25)
    fig.suptitle('VZ = WZ + ZZ; error bands are statistical only', fontsize=9)
    save_figure(fig, output, f'{cat}_{tag}_{ch}_{target}_mass')


def yield_plots(rows, output):
    for ch in CHANNELS+['combined']:
        for target in ['SR', 'VR']:
            subset = [r for r in rows if r['channel'] == ch and r['target'] == target]
            x = np.arange(len(subset))
            labels = [r['category']+'\n'+r['tag'] for r in subset]
            fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
            styles = [
                ('obs_nozh', 'black', 'o', 'Target Data − TT/WW/tZq − VZ'),
                ('obs_withzh', 'gray', 'o', 'Target additionally subtracting ZHbb'),
                ('pred_nozh', 'tab:blue', 's', 'Alpha estimate without ZH subtraction'),
                ('pred_withzh', 'tab:green', '^', 'Alpha estimate with ZH subtraction'),
                ('pred_total_withzh', 'tab:red', 'd', 'Corrected DY + target ZHbb')]
            for i, (key, color, marker, label) in enumerate(styles):
                ax.errorbar(x+(i-2)*.12, [r[key] for r in subset],
                    [r[key+'_stat'] for r in subset], fmt=marker, color=color, ms=5, label=label)
            ax.set_yscale('symlog', linthresh=10)
            ax.set_xticks(x, labels); ax.set_ylabel('Yield, m(2ℓ2q) ≥300 GeV')
            ax.set_title(f'2018 {ch}, {target} — all standard categories; statistical errors only')
            ax.legend(fontsize=8); ax.grid(axis='y', alpha=.2)
            save_figure(fig, output, f'{ch}_{target}_yield_summary')
            fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
            for key, label, color, offset in [
                ('fraction_dy_prediction_removed', 'ZH transferred from CR / old DY estimate', 'tab:green', -.16),
                ('fraction_target_zh', 'Target ZH / target Data−other backgrounds', 'tab:purple', 0),
                ('fraction_total_prediction_change', 'Corrected DY+ZH change / old estimate', 'tab:red', .16)]:
                ax.plot(x+offset, [100*r[key] for r in subset], 'o', color=color, label=label)
            ax.axhline(0, color='gray', ls='--'); ax.set_xticks(x, labels)
            ax.set_ylabel('Fraction [%]'); ax.set_title(f'2018 {ch}, {target} — ZHbb impact (central values)')
            ax.legend(fontsize=8); ax.grid(axis='y', alpha=.2)
            save_figure(fig, output, f'{ch}_{target}_zh_impact_summary')


def hadronic_mass_plots(arrays, cat, tag, output, records):
    """Hadronic shape diagnostic only: alpha does not predict target q shape."""
    qedges = np.array([40, 55, 70, 85, 105, 120, 135, 150, 165, 180.])
    per_channel = []
    for ch in CHANNELS:
        h = {}
        for proc in PROCESSES:
            h[proc] = np.zeros((2, len(qedges)-1))
            for reg in ['CR', 'SR', 'VR']:
                a = arrays[ch, reg, proc]
                h[proc] += hist({'mass': a['mass2q'], 'weight': a['weight']}, qedges)
        per_channel.append(h)
    combined = {proc: sum(h[proc] for h in per_channel) for proc in PROCESSES}
    for ch, h in list(zip(CHANNELS, per_channel))+[('combined', combined)]:
        t = h['Data'][0]-h['TTplusWW'][0]-h['VZ'][0]
        vt = h['Data'][1]+h['TTplusWW'][1]+h['VZ'][1]
        zh, vzh = h['ZHbb']; dy, vdy = h['DY']
        x = .5*(qedges[:-1]+qedges[1:]); widths = np.diff(qedges)
        fig, (ax, ar) = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                            gridspec_kw={'height_ratios': [3, 1.4]}, constrained_layout=True)
        ax.errorbar(x-0.5, t/widths, np.sqrt(vt)/widths, fmt='o', color='black', ms=4,
                    label='Data − TT/WW/tZq − VZ')
        ax.errorbar(x+0.5, (t-zh)/widths, np.sqrt(vt+vzh)/widths, fmt='o',
                    color='tab:green', mfc='white', ms=4, label='Additionally subtract ZHbb')
        ax.stairs(dy/widths, qedges, color='tab:orange', label='Direct DY MC')
        ax.stairs(zh/widths, qedges, color='tab:purple', label='ZHbb MC')
        ax.stairs((dy+zh)/widths, qedges, color='tab:red', ls='--', label='DY + ZHbb MC')
        ax.set_ylabel('Events / GeV'); ax.legend(fontsize=8)
        ax.set_title(f'2018 {cat}, {tag}, {ch} — m(2ℓ2q) ≥300 GeV\nHadronic-mass diagnostic; not an alpha prediction of the q-mass shape')
        for num, var, color, label, off in [(t, vt, 'black', '(Data−other BG) / DY', -.5),
                (t-zh, vt+vzh, 'tab:green', '(Data−other BG−ZHbb) / DY', .5)]:
            ar.errorbar(x+off, divide(num, dy), divide(np.sqrt(var), dy), fmt='o', color=color, ms=4, label=label)
        band = divide(np.sqrt(vdy), dy)
        ar.fill_between(qedges, np.r_[1-band, (1-band)[-1]], np.r_[1+band, (1+band)[-1]],
                        step='post', color='tab:orange', alpha=.15)
        ar.axhline(1, color='gray', ls='--'); ar.set_ylim(-.5, 2.5)
        ar.set_ylabel('Signed data / DY'); ar.set_xlabel('Hadronic mass m(qq) or m(J) [GeV]')
        ar.legend(fontsize=8)
        for axis in [ax, ar]:
            axis.axvspan(70, 105, color='tab:green', alpha=.05)
            axis.axvspan(105, 135, color='tab:red', alpha=.07)
            axis.axvspan(135, 165, color='tab:purple', alpha=.05)
        save_figure(fig, output, f'{cat}_{tag}_{ch}_mass2q')
        for i in range(len(qedges)-1):
            row = dict(category=cat, tag=tag, channel=ch, low=qedges[i], high=qedges[i+1])
            for proc in PROCESSES:
                row[proc] = h[proc][0, i]; row[proc+'_variance'] = h[proc][1, i]
            records.append(row)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base', type=Path, required=True)
    p.add_argument('--cr-base', type=Path, required=True)
    p.add_argument('--zh-base', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--fine-binning', action='store_true', help='DY-only finer transfer-binning sensitivity check')
    p.add_argument('--no-plots', action='store_true', help='Save numerical results only')
    args = p.parse_args()
    if args.fine_binning:
        INITIAL_EDGES['resolved'] = list(np.arange(300., 1101., 50.))+[1400., 1800., 2300., 3500., np.inf]
        INITIAL_EDGES['merged'] = list(np.arange(300., 1101., 100.))+[1400., 1800., 3500., np.inf]
    if args.output.exists(): raise FileExistsError('Use a new study-output directory')
    args.output.mkdir(parents=True)
    audit, rows, bins, qbins, edge_info, archive = [], [], [], [], [], {}
    for cat in CATEGORIES:
        for tag in TAGS:
            print('Comparing', cat, tag, flush=True)
            arrays = read_inputs(args, cat, tag, audit)
            if not args.no_plots: hadronic_mass_plots(arrays, cat, tag, args.output, qbins)
            edges, merges = adaptive_edges(arrays, cat)
            edge_info.append(dict(category=cat, tag=tag, edges=[None if np.isinf(e) else float(e) for e in edges],
                                  removed_edges=merges, criterion='Positive DY sumw in CR/SR/VR separately in ee and mumu; never uses data or ZH'))
            for target in ['SR', 'VR']:
                results = [derive(make_histograms(arrays, ch, target, edges)) for ch in CHANNELS]
                for ch, r in list(zip(CHANNELS, results))+[('combined', combine(results))]:
                    row = yield_row(r, cat, tag, ch, target); rows.append(row)
                    for i in range(len(edges)-1):
                        b = dict(category=cat, tag=tag, channel=ch, target=target,
                                 low=float(edges[i]), high=float(edges[i+1]))
                        for key, value in r.items():
                            name = '_'.join(key) if isinstance(key, tuple) else key
                            if isinstance(key, str) and key.startswith('cov_'):
                                b[name] = float(value[i]); continue
                            b[name] = float(value[0, i]); b[name+'_variance'] = float(value[1, i])
                        bins.append(b)
                    prefix = '__'.join([cat, tag, ch, target])
                    archive[prefix+'__edges'] = edges
                    for key, value in r.items():
                        archive[prefix+'__'+('_'.join(key) if isinstance(key, tuple) else key)] = value
                    if not args.no_plots: distribution_plot(r, edges, cat, tag, ch, target, args.output)
                    print({k: row[k] for k in ['category','tag','channel','target','obs_nozh','pred_nozh','pred_withzh','target_ZHbb','transferred_zh','closure_nozh','closure_withzh']}, flush=True)
            del arrays
    if not args.no_plots: yield_plots(rows, args.output)
    for filename, records in [('yields.csv', rows), ('bins.csv', bins), ('mass2q_bins.csv', qbins)]:
        if not records: continue
        with (args.output/filename).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)
    np.savez_compressed(args.output/'histograms.npz', **archive)
    (args.output/'input_audit.json').write_text(json.dumps(audit, indent=2))
    (args.output/'binning.json').write_text(json.dumps(edge_info, indent=2))
    metadata = dict(alpha='Nominal DY target/CR, recomputed in documented coarse positive-yield bins',
        mass='All m(2l2q)>=300 including infinite high-mass tail; no control region in this mass',
        cr='Hadronic mass 40–180 outside SR 70–105 and VR 135–165, retaining 105–135',
        backgrounds='TTplusWW is TTTo2L2Nu+WWTo2L2Nu+tZq; VZ is WZTo2Q2L+ZZTo2Q2L',
        before='P0=alpha*(Data_CR-TTplusWW_CR-VZ_CR); O0=Data_target-TTplusWW_target-VZ_target',
        after='P1=alpha*(Data_CR-TTplusWW_CR-VZ_CR-ZHbb_CR); O1=O0-ZHbb_target',
        fair_total_comparison='Compare P1+ZHbb_target and P0 to the SAME O0 target curve',
        combined='Sum channel-specific predictions/variances; combined DY ratio saved for reference, not applied',
        uncertainties='Statistical sumw2 and finite independent DY numerator/denominator; luminosity, SF, normalization and modelling systematics omitted',
        correlation='Before/after are correlated; delta variance uses alpha*ZH_CR and target ZH, not quadrature of two full predictions',
        chi2='Approximate Gaussian diagnostic only; sparse/negative signed bins and systematic correlations require a proper likelihood',
        only_ZHbb='Does not include ggZH, ZHcc or other Higgs processes',
        note='Baseline cached reduced MC weights are the existing analysis convention, not independently recalibrated in this study',
        fine_binning=args.fine_binning,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), plots=0 if args.no_plots else 2*66, rows=len(rows))
    (args.output/'metadata.json').write_text(json.dumps(metadata, indent=2))
    print('Study complete:', args.output, flush=True)


if __name__ == '__main__': main()
