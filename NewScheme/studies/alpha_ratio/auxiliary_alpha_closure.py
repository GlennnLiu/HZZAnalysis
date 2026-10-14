"""Orthogonal MET-selection alpha closure, using original skims (not RedTrees).

This is a diagnostic with freshly computed, coarse-bin nominal DY ratios,
not an application or replacement of the production alpha histograms.
All m(2l2q) >= 300 GeV is retained, including a >=3500 overflow bin.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
from concurrent.futures import ThreadPoolExecutor

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'
import awkward as ak
import numpy as np
import uproot
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

Q = np.array([40, 55, 70, 85, 105, 120, 135, 150, 165, 180.])
X = {'resolved': np.array([300, 400, 500, 650, 850, 1100, 1600, 3500, np.inf]),
     'merged': np.array([300, 500, 850, 3500, np.inf])}
MET = np.array([50, 100, 150, 250.])
SELECTIONS = {'nominal': [0, 1, 2], 'MET_0_50': [0], 'MET_50_100': [1],
              'MET_100_150': [2], 'MET_150_250': [3], 'MET_gt250': [4],
              'MET_gt150': [3, 4]}
REGIONS = {'CR_broad_exclVR': [0, 1, 4, 5, 8],
           'CR_broad_inclVR': [0, 1, 4, 5, 6, 7, 8],
           'CR_noH_exclVR': [0, 1, 8], 'SR': [2, 3], 'VR': [6, 7]}
PROCESSES = ['Data', 'DY', 'TTplusWW', 'VZ']


class AuxiliaryDeepJet:
    """Reproduce btageff.py weights in memory, validating against nominal entries.

    Includes its legacy second-jet charm lookup in the both-fail case, solely
    for consistency with stored weights; this is not a new SF prescription.
    Never changes input trees or calibration files.
    """
    def __init__(self, args, sample):
        import correctionlib
        self.corrections = correctionlib.CorrectionSet.from_file(
            str(args.deepjet_dir / 'btagging.json.gz'))
        eff_sample = 'TTTo2L2Nu' if sample == 'WWTo2L2Nu' else sample
        self.eff_path = args.deepjet_dir / ('btag_eff_%s.root' % eff_sample)
        with uproot.open(self.eff_path) as f:
            self.eff = {flavor: (f[key].values(flow=True),
                                f[key].axes[0].edges(), f[key].axes[1].edges())
                        for flavor, key in [(0, 'eff_udsg'), (4, 'eff_c'), (5, 'eff_b')]}
        self.validation_max = 0.; self.validated = 0; self.repaired = 0

    def calculate(self, a, mask):
        jets, flavors = [], []
        both_fail = (a['jet_1_deepbtag'][mask] < .2783) & (a['jet_2_deepbtag'][mask] < .2783)
        for j in [1, 2]:
            flavor = a['jet_%d_hadronflavor' % j][mask].astype(np.int64)
            pt = a['jet_%d_pt' % j][mask].astype(float)
            eta = np.abs(a['jet_%d_eta' % j][mask].astype(float))
            sf = np.empty(len(pt)); efficiency = np.empty(len(pt))
            if not np.all(np.isin(flavor, [0, 4, 5])):
                raise ValueError('Unexpected jet hadron flavour')
            eff_flavor = np.where(both_fail & (flavor == 4), 5, flavor) if j == 2 else flavor
            for flav in [0, 4, 5]:
                m = flavor == flav
                if m.any():
                    correction = self.corrections['deepJet_incl' if flav == 0 else 'deepJet_comb']
                    sf[m] = correction.evaluate('central', 'M', flavor[m], eta[m], pt[m])
                m = eff_flavor == flav
                if m.any():
                    values, xe, ye = self.eff[flav]
                    efficiency[m] = values[np.searchsorted(xe, pt[m], side='right'),
                                           np.searchsorted(ye, eta[m], side='right')]
            passed = a['jet_%d_deepbtag' % j][mask] >= .2783
            if np.any((~passed) & (efficiency >= 1)):
                raise ValueError('Invalid fail-jet efficiency')
            per_jet = sf.copy()
            per_jet[~passed] = (1-sf[~passed]*efficiency[~passed])/(1-efficiency[~passed])
            jets.append(per_jet)
        return jets[0]*jets[1]

    def extend(self, a, selected):
        out = a['w_deepjet'].astype(float).copy()
        if not selected.any():
            return out
        computed = self.calculate(a, selected)
        stored = out[selected]
        nominal = a['Met'][selected] <= 150
        if np.any(nominal & (stored < -900)):
            raise ValueError('Nominal selected DeepJet sentinel')
        if nominal.any():
            delta = np.abs(computed[nominal]-stored[nominal])/np.maximum(np.abs(stored[nominal]), 1e-12)
            self.validation_max = max(self.validation_max, float(delta.max()))
            self.validated += int(nominal.sum())
            if self.validation_max > 2e-5:
                raise ValueError('Recomputed SF does not reproduce stored nominal weights: %g' % self.validation_max)
        missing = stored < -900
        if np.any(missing & nominal):
            raise ValueError('Refusing to replace nominal weights')
        indices = np.flatnonzero(selected)[missing]
        out[indices] = computed[missing]
        self.repaired += len(indices)
        return out


def scan(args, repair=False):
    sys.path.insert(0, str(args.analysis / 'FitTemplate'))
    import make_reduced_trees as producer
    config_path = args.analysis / 'cards' / ('config_UL%s.yml' % args.era)
    config = yaml.safe_load(config_path.read_text())
    sf_path = producer.DEFAULT_SF_DIR / ('NetSF_signal_%s.yml' % args.era)
    sfs = yaml.safe_load(sf_path.read_text())
    hist = {cat: {proc: np.zeros((2, 5, len(X[cat])-1, len(Q)-1, 2))
                  for proc in PROCESSES} for cat in X}
    previous_audit = {}
    if repair:
        with np.load(args.output / 'source_histograms.npz') as f:
            for cat in X:
                for proc in ['DY', 'Data']:
                    hist[cat][proc] = f['%s__%s' % (cat, proc)]
        previous_audit = {r['sample']: r for r in json.loads(
            (args.output / 'source_audit.json').read_text())['samples']}
    audit = []
    executor = ThreadPoolExecutor(max_workers=2)
    for sample in config['Samples_lists']:
        group = producer.sample_group(sample)
        if repair and group in ['DY', 'Data']:
            audit.append(previous_audit[sample]); continue
        path = producer.input_path(config, args.era, sample)
        needed = {'Met', 'foundZ1LCandidate', 'foundTTCRCandidate', 'isEE', 'isMuMu',
                  'foundZ2JCandidate', 'foundZ2MergedCandidata', 'mass2l2jet',
                  'mass2jet', 'mass2lj', 'massmerged', 'KD_jjVBF', 'KD_jVBF',
                  'KD_JVBF', 'KD_JJVBF', 'jet_1_deepbtag', 'jet_2_deepbtag',
                  'particleNetZvsQCD', 'particleNetZbbvslight'}
        if group != 'Data':
            needed.add('EventWeight')
            if group != 'DY':
                needed.add('w_deepjet')
                needed.update({'jet_%d_%s' % (j, variable) for j in [1, 2]
                               for variable in ['pt', 'eta', 'hadronflavor']})
            if sample in producer.PARTICLENET_SAMPLES:
                needed.update({'isbjet', 'iscjet', 'islightjet', 'ptmerged'})
        row = dict(sample=sample, path=str(path), size=path.stat().st_size,
                   mtime=path.stat().st_mtime, entries=0, falseZ=0,
                   falseZ_sameflavour=0, falseZ_topflag=0)
        print('Scanning %s' % sample, flush=True)
        deepjet = (AuxiliaryDeepJet(args, sample)
                   if group not in ['Data', 'DY'] and not args.auxiliary_unit_deepjet else None)
        with uproot.open(path, decompression_executor=executor,
                         interpretation_executor=uproot.TrivialExecutor()) as f:
            tree = f['passedEvents']
            missing = needed - set(tree.keys())
            if missing:
                raise KeyError('%s: missing %s' % (path, sorted(missing)))
            norm = None if group == 'Data' else float(f['sumWeights'].values().sum())
            if group != 'Data' and (not np.isfinite(norm) or norm == 0):
                raise ValueError('Invalid sumWeights for %s' % sample)
            row['sumWeights'] = norm
            for events in tree.iterate(sorted(needed), step_size=args.chunk_size,
                                       library='ak'):
                a = {key: ak.to_numpy(events[key]) for key in needed}
                row['entries'] += len(events)
                falsez = ~a['foundZ1LCandidate'].astype(bool)
                row['falseZ'] += int(falsez.sum())
                row['falseZ_sameflavour'] += int((falsez & (a['isEE'] | a['isMuMu'])).sum())
                row['falseZ_topflag'] += int((falsez & a['foundTTCRCandidate']).sum())
                # side='left': nominal includes MET=150, auxiliary is strictly >150.
                mi = np.searchsorted(MET, a['Met'], side='left')
                for cat in X:
                    if cat == 'resolved':
                        mass, q = a['mass2l2jet'], a['mass2jet']
                        base = (a['foundZ2JCandidate'] & (a['KD_jjVBF'] < .5)
                                & (a['KD_jVBF'] < .7)
                                & ((a['jet_1_deepbtag'] < .2783)
                                   | (a['jet_2_deepbtag'] < .2783)))
                    else:
                        mass, q = a['mass2lj'], a['massmerged']
                        base = (a['foundZ2MergedCandidata'] & (a['KD_JJVBF'] < .5)
                                & (a['KD_JVBF'] < .7) & (a['particleNetZvsQCD'] > .9)
                                & (a['particleNetZbbvslight'] < .9))
                    nx, nq = len(X[cat])-1, len(Q)-1
                    xi = np.searchsorted(X[cat], mass, side='right') - 1
                    qi = np.searchsorted(Q, q, side='right') - 1
                    base &= (a['foundZ1LCandidate'] & np.isfinite(mass)
                             & np.isfinite(q) & np.isfinite(a['Met'])
                             & (xi >= 0) & (xi < nx) & (qi >= 0) & (qi < nq))
                    weights = (np.ones(len(events)) if group == 'Data' else
                               producer.mc_weights(events, config, args.era, sample,
                                                   group, cat, 'untag', norm, sfs))
                    if cat == 'resolved' and deepjet is not None:
                        selected = base & (a['isEE'] | a['isMuMu'])
                        dj = deepjet.extend(a, selected)
                        weights = (a['EventWeight'] * config['lumi'][args.era] * 1000
                                   * config['samples_inf'][sample][1] / norm * dj)
                    elif cat == 'resolved' and group not in ['Data', 'DY']:
                        # Explicit feasibility approximation, never mistaken for calibrated output.
                        selected = base & (a['isEE'] | a['isMuMu'])
                        if np.any(selected & (a['Met'] <= 150) & (a['w_deepjet'] < -900)):
                            raise ValueError('Invalid nominal DeepJet weight')
                        dj = np.where(a['Met'] > 150, 1., a['w_deepjet'])
                        weights = (a['EventWeight'] * config['lumi'][args.era] * 1000
                                   * config['samples_inf'][sample][1] / norm * dj)
                    if not np.all(np.isfinite(weights[base])):
                        raise ValueError('Nonfinite selected weights: %s %s' % (sample, cat))
                    for channel, field in enumerate(['isEE', 'isMuMu']):
                        mask = base & a[field]
                        indices = mi[mask]*nx*nq + xi[mask]*nq + qi[mask]
                        for power in [1, 2]:
                            hist[cat][group][channel, ..., power-1] += np.bincount(
                                indices, weights=weights[mask]**power,
                                minlength=5*nx*nq).reshape(5, nx, nq)
        if deepjet is not None:
            row['deepjet_extension'] = dict(repaired=deepjet.repaired,
                                            validated=deepjet.validated,
                                            nominal_max_relative_difference=deepjet.validation_max,
                                            efficiency_file=str(deepjet.eff_path))
            print('  DeepJet validated=%d repaired=%d max_relative_difference=%g' %
                  (deepjet.validated, deepjet.repaired, deepjet.validation_max), flush=True)
        elif group not in ['Data', 'DY']:
            row['deepjet_extension'] = 'High-MET SF=1: uncalibrated feasibility study only'
        audit.append(row)
        print('  entries=%d; falseZ=%d; falseZ same-flavour=%d' %
              (row['entries'], row['falseZ'], row['falseZ_sameflavour']), flush=True)
    arrays = {'%s__%s' % (cat, proc): values
              for cat, components in hist.items() for proc, values in components.items()}
    np.savez_compressed(args.output / 'source_histograms.npz', **arrays)
    (args.output / 'source_audit.json').write_text(json.dumps(
        dict(era=args.era, config=config, sf_file=str(sf_path), scale_factors=sfs,
             samples=audit, mass_edges={k: [str(v) for v in values] for k, values in X.items()},
             q_edges=Q.tolist(), MET_upper_edges=MET.tolist(),
             deepjet_auxiliary_extension=('High-MET SF=1: uncalibrated feasibility study only'
                                         if args.auxiliary_unit_deepjet else
                                         'Validated against stored nominal weights; legacy btageff.py formula'),
             deepjet_calibration=str(args.deepjet_dir / 'btagging.json.gz')), indent=2))
    executor.shutdown()
    return hist


def region(components, indices, tt_scale=1., vz_scale=1.):
    p = {name: values[:, indices, :].sum(axis=1) for name, values in components.items()}
    z = p['Data'][:, 0] - tt_scale*p['TTplusWW'][:, 0] - vz_scale*p['VZ'][:, 0]
    variance = (p['Data'][:, 1] + tt_scale**2*p['TTplusWW'][:, 1]
                + vz_scale**2*p['VZ'][:, 1])
    return z, variance, p['DY'][:, 0], p['DY'][:, 1], p


def closure(c, t):
    zc, vc, dc, vdc = c[:4]
    zt, vt, dt, vdt = t[:4]
    with np.errstate(divide='ignore', invalid='ignore'):
        alpha = dt / dc
        pred = alpha * zc
        predvar = alpha**2*vc + (zc/dc)**2*vdt + (zc*dt/dc**2)**2*vdc
        k = zt / pred
        kv = vt/pred**2 + zt**2*predvar/pred**4
        valid = (zc > 0) & (dc > 0) & (dt > 0) & (pred > 0)
        return [np.where(valid, v, np.nan) for v in (alpha, pred, predvar, k, kv)]


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def save(fig, path):
    for suffix in ['pdf', 'png']:
        fig.savefig(str(path) + '.' + suffix, dpi=170, bbox_inches='tight')
    plt.close(fig)


def report(args, hist):
    metadata = json.loads((args.output / 'source_audit.json').read_text())
    provisional = 'SF=1' in metadata['deepjet_auxiliary_extension']
    weight_note = ('\nHigh-MET non-DY DeepJet SF=1: feasibility approximation, NOT calibrated'
                   if provisional else '')
    yields, ratios, stresses = [], [], []
    for cat in X:
        for channel in ['ee', 'mumu', 'combined']:
            components = {proc: (values.sum(axis=0) if channel == 'combined'
                                 else values[['ee', 'mumu'].index(channel)])
                          for proc, values in hist[cat].items()}
            selections = {label: {p: v[metbins].sum(axis=0) for p, v in components.items()}
                          for label, metbins in SELECTIONS.items()}
            measures = {}
            for selection, cells in selections.items():
                rs = {name: region(cells, indices) for name, indices in REGIONS.items()}
                for rname, r in rs.items():
                    z, vz, d, vd, p = r
                    for i in range(len(z)+1):
                        def val(a): return float(a.sum() if i == len(z) else a[i])
                        yields.append(dict(category=cat, channel=channel, selection=selection,
                                           region=rname, xlow=300 if i == len(z) else X[cat][i],
                                           xhigh=np.inf if i == len(z) else X[cat][i+1],
                                           integrated=i == len(z), signed=val(z),
                                           signed_stat=np.sqrt(val(vz)), DY=val(d),
                                           DY_stat=np.sqrt(val(vd)), Data=val(p['Data'][:, 0]),
                                           TTplusWW=val(p['TTplusWW'][:, 0]), VZ=val(p['VZ'][:, 0])))
                for cname in ['CR_broad_exclVR', 'CR_broad_inclVR', 'CR_noH_exclVR']:
                    for target in ['SR', 'VR']:
                        if target == 'VR' and cname == 'CR_broad_inclVR':
                            continue  # avoids target/source overlap and wrong covariance
                        c, t = rs[cname], rs[target]
                        result = closure(c, t)
                        measures[(selection, cname, target)] = result
                        for integrated in [False, True]:
                            cc = [np.array([a.sum()]) for a in c[:4]] if integrated else c
                            tt = [np.array([a.sum()]) for a in t[:4]] if integrated else t
                            aa, pp, vp, kk, vk = closure(cc, tt)
                            for i in range(len(kk)):
                                eff = lambda a, v: float(a*a/v) if v > 0 else 0.
                                ratios.append(dict(category=cat, channel=channel,
                                                   selection=selection, source=cname, target=target,
                                                   integrated=integrated,
                                                   xlow=300 if integrated else X[cat][i],
                                                   xhigh=np.inf if integrated else X[cat][i+1],
                                                   alpha=aa[i], prediction=pp[i],
                                                   prediction_stat=np.sqrt(vp[i]), kappa=kk[i],
                                                   kappa_stat=np.sqrt(vk[i]),
                                                   source_signed_neff=eff(cc[0][i], cc[1][i]),
                                                   target_signed_neff=eff(tt[0][i], tt[1][i]),
                                                   source_DY_neff=eff(cc[2][i], cc[3][i]),
                                                   target_DY_neff=eff(tt[2][i], tt[3][i])))
                for tt_scale, vz_scale in [(0.8, 1.), (1.2, 1.), (1., .7), (1., 1.3)]:
                    c = region(cells, REGIONS['CR_broad_exclVR'], tt_scale, vz_scale)
                    t = region(cells, REGIONS['SR'], tt_scale, vz_scale)
                    result = closure([np.array([a.sum()]) for a in c[:4]],
                                     [np.array([a.sum()]) for a in t[:4]])
                    stresses.append(dict(category=cat, channel=channel, selection=selection,
                                         TT_scale=tt_scale, VZ_scale=vz_scale,
                                         kappa=result[3][0], kappa_stat=np.sqrt(result[4][0])))
            cname = 'CR_broad_exclVR'
            bins = np.arange(len(X[cat])-1)
            labels = ['%g-%g' % (a, b) if np.isfinite(b) else '>=%g' % a
                      for a, b in zip(X[cat][:-1], X[cat][1:])]
            for target in ['SR', 'VR']:
                fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
                for offset, selection in [(-.15, 'nominal'), (0., 'MET_150_250'), (.15, 'MET_gt250')]:
                    r = measures[(selection, cname, target)]
                    axes[0].errorbar(bins+offset, r[3], np.sqrt(r[4]), fmt='o', label=selection)
                nom, aux = [measures[(s, cname, target)] for s in ['nominal', 'MET_gt150']]
                with np.errstate(divide='ignore', invalid='ignore'):
                    transfer = nom[3]/aux[3]
                    error = np.sqrt(nom[4]/aux[3]**2 + nom[3]**2*aux[4]/aux[3]**4)
                axes[1].errorbar(bins, transfer, error, fmt='o', color='black')
                for ax in axes:
                    ax.axhline(1, color='gray', ls='--'); ax.grid(alpha=.15)
                axes[0].set_ylabel('kappa = signed target / alpha prediction')
                axes[0].set_ylim(-1, 3); axes[0].legend()
                axes[1].set_ylabel('kappa(nominal) / kappa(MET > 150)')
                axes[1].set_ylim(-1, 3); axes[1].set_xticks(bins)
                axes[1].set_xticklabels(labels, rotation=30)
                axes[1].set_xlabel('m(2l2q) [GeV]; all bins are search-domain bins')
                fig.suptitle('%s %s untag, %s, target %s\nStatistical errors only; CR=40-70,105-135,165-180 GeV' %
                             (args.era, cat, channel, target) + (weight_note if cat == 'resolved' else ''))
                save(fig, args.output / ('%s_%s_%s_MET_transfer' % (cat, channel, target)))
            fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
            nom, aux = [measures[(s, cname, 'SR')] for s in ['nominal', 'MET_gt150']]
            for ax, selection, result in zip(axes[:2], ['nominal', 'MET_gt150'], [nom, aux]):
                target = region(selections[selection], REGIONS['SR'])
                ax.errorbar(bins-.06, target[0], np.sqrt(target[1]), fmt='o', label='Data - TTplusWW - VZ')
                ax.errorbar(bins+.06, result[1], np.sqrt(result[2]), fmt='s', label='DY alpha x signed CR')
                ax.set_ylabel('Yield / mass bin'); ax.set_title(selection); ax.legend()
            target = region(selections['nominal'], REGIONS['SR'])
            corrected = aux[3]*nom[1]
            variance = aux[3]**2*nom[2] + nom[1]**2*aux[4]
            axes[2].errorbar(bins-.06, target[0], np.sqrt(target[1]), fmt='o', label='Nominal signed SR (diagnostic only)')
            axes[2].errorbar(bins+.06, corrected, np.sqrt(variance), fmt='s', label='Prediction x kappa(MET > 150)')
            axes[2].set_ylabel('Yield / mass bin'); axes[2].legend()
            axes[2].set_xticks(bins); axes[2].set_xticklabels(labels, rotation=30)
            axes[2].set_xlabel('m(2l2q) [GeV]')
            fig.suptitle('%s %s untag, %s: independent-selection SR closure\nFresh coarse-bin alpha; statistical errors only; no production changes' %
                         (args.era, cat, channel) + (weight_note if cat == 'resolved' else ''))
            save(fig, args.output / ('%s_%s_SR_auxiliary_closure' % (cat, channel)))
    write_csv(args.output / 'yields.csv', yields)
    write_csv(args.output / 'alpha_closure.csv', ratios)
    write_csv(args.output / 'subtraction_stress.csv', stresses)
    print('Saved auxiliary-selection study to %s' % args.output, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis', type=Path, default=Path(__file__).resolve().parents[3])
    p.add_argument('--era', choices=['2018'], default='2018')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--chunk-size', type=int, default=200000)
    p.add_argument('--reuse-cache', action='store_true', help='Explicitly reuse saved histograms; does not recheck source freshness.')
    p.add_argument('--repair-auxiliary-deepjet', action='store_true', help='Rescan non-DY MC, retaining cached DY and Data.')
    p.add_argument('--auxiliary-unit-deepjet', action='store_true', help='Explicit uncalibrated feasibility mode: resolved high-MET non-DY DeepJet SF=1.')
    p.add_argument('--deepjet-dir', type=Path, default=Path('/cms/user/guojl/ME_test/CMSSW_10_6_26/src/HZZAnalysis/deepjetSFsfiles/2018'))
    args = p.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    if args.reuse_cache:
        metadata = json.loads((args.output / 'source_audit.json').read_text())
        if 'deepjet_auxiliary_extension' not in metadata:
            raise ValueError('Cache has unvalidated high-MET DeepJet weights; run --repair-auxiliary-deepjet')
        with np.load(args.output / 'source_histograms.npz') as f:
            hist = {cat: {proc: f['%s__%s' % (cat, proc)] for proc in PROCESSES} for cat in X}
    else:
        hist = scan(args, repair=args.repair_auxiliary_deepjet)
    report(args, hist)


if __name__ == '__main__':
    main()
