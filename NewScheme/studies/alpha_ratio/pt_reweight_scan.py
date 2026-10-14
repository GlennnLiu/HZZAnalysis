"""Cache nominal 2018 (search mass, hadronic pT) sumw/sumw2, without ZH."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import re
import sys

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import awkward as ak
import numpy as np
import uproot
import yaml

X = {'resolved': np.array([300,400,500,600,750,900,1100,1400,1800,2300,3500,np.inf]),
     'merged': np.array([300,500,750,1100,1600,3500,np.inf])}
PT = np.array([0,30,60,100,150,200,250,300,350,400,500,650,850,1200,np.inf])
REGIONS = ('CR', 'SR', 'VR')
CHANNELS = ('ee', 'mumu')
TAGS = ('btag', 'untag', 'vbftag')
PROCESSES = ('Data', 'DY', 'TTplusWW', 'VZ')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis', type=Path, default=Path('/cms/user/liugeliang/CMS/HZZAnalysis'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--chunk-size', default='60 MB')
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new cache directory')
    sys.path.insert(0, str(args.analysis/'FitTemplate'))
    import make_reduced_trees as producer
    config_path = args.analysis/'cards/config_UL2018.yml'
    config = yaml.safe_load(config_path.read_text())
    sf_path = producer.DEFAULT_SF_DIR/'NetSF_signal_2018.yml'
    sfs = yaml.safe_load(sf_path.read_text())
    expressions, needed = {}, {'pt2jet', 'ptmerged'}
    for cat in X:
        needed.update(producer.BRANCHES[cat].values())
        needed.update(producer.selection_branches(config, ['SR'], CHANNELS, cat, TAGS))
        qname = producer.BRANCHES[cat]['mass2q']
        for ch in CHANNELS:
            for tag in TAGS:
                expression = config['cut']['SR'][producer.CHANNELS[ch]][cat][tag]
                expression, count = re.subn(r'\(\('+qname+r'>70\)\s*&\s*\('+qname+r'<105\)\)',
                    '(('+qname+'>40) & ('+qname+'<180))', expression)
                if count != 1:
                    raise ValueError('Unrecognized hadronic-mass selection: '+expression)
                expressions[cat, tag, ch] = expression
    cache = {f'{cat}__{tag}__{proc}': np.zeros((2,3,len(X[cat])-1,len(PT)-1,2))
             for cat in X for tag in TAGS for proc in PROCESSES}
    audit = []
    args.output.mkdir(parents=True)
    with ThreadPoolExecutor(max_workers=2) as executor:
        for sample in config['Samples_lists']:
            group = producer.sample_group(sample)
            path = producer.input_path(config, '2018', sample)
            fields = set(needed)
            if group != 'Data': fields.add('EventWeight')
            if group not in ('Data', 'DY'): fields.add('w_deepjet')
            if sample in producer.PARTICLENET_SAMPLES:
                fields.update(['isbjet','iscjet','islightjet','ptmerged'])
            row = dict(sample=sample, group=group, path=str(path), bytes=path.stat().st_size,
                       mtime=path.stat().st_mtime, entries_read=0, selected={})
            print('Scanning '+sample, flush=True)
            with uproot.open(path, decompression_executor=executor,
                             interpretation_executor=uproot.TrivialExecutor()) as f:
                tree = f['passedEvents']
                missing = fields-set(tree.keys())
                if missing: raise KeyError(f'{sample}: {sorted(missing)}')
                norm = None if group=='Data' else float(f['sumWeights'].values().sum())
                if group!='Data' and (not np.isfinite(norm) or norm==0):
                    raise ValueError('Invalid MC normalization')
                row['normalization'] = norm
                for events in tree.iterate(sorted(fields), step_size=args.chunk_size, library='ak'):
                    row['entries_read'] += len(events)
                    for cat in X:
                        mass = ak.to_numpy(events[producer.BRANCHES[cat]['mass']])
                        q = ak.to_numpy(events[producer.BRANCHES[cat]['mass2q']])
                        pt = ak.to_numpy(events['pt2jet' if cat=='resolved' else 'ptmerged'])
                        sr = (q>70)&(q<105)
                        vr = (q>135)&(q<165)
                        cr = ((q<70)|(q>105))&~vr
                        xi = np.searchsorted(X[cat], mass, side='right')-1
                        pi = np.searchsorted(PT, pt, side='right')-1
                        nx, np_ = len(X[cat])-1, len(PT)-1
                        for tag in TAGS:
                            for ci, ch in enumerate(CHANNELS):
                                selected = ak.to_numpy(ak.numexpr.evaluate(expressions[cat,tag,ch], events))
                                selected &= mass>=300
                                if np.any(selected & (~np.isfinite(mass)|~np.isfinite(q)|~np.isfinite(pt)|(pt<0))):
                                    raise ValueError('Invalid selected kinematics')
                                chosen = events[selected]
                                if cat=='resolved' and tag!='vbftag' and group not in ('Data','DY'):
                                    if np.any(ak.to_numpy(chosen['w_deepjet']) < -900):
                                        raise ValueError('DeepJet sentinel in nominal selection')
                                weight = (np.ones(len(chosen)) if group=='Data' else
                                          producer.mc_weights(chosen, config, '2018', sample, group, cat, tag, norm, sfs))
                                if not np.all(np.isfinite(weight)):
                                    raise ValueError('Nonfinite selected weights')
                                if np.any(xi[selected]<0)|np.any(pi[selected]<0):
                                    raise ValueError('Selected event outside histogram domain')
                                for ri, regionmask in enumerate((cr, sr, vr)):
                                    part = regionmask[selected]
                                    indices = xi[selected][part]*np_+pi[selected][part]
                                    values = cache[f'{cat}__{tag}__{group}'][ci,ri]
                                    for power in (1,2):
                                        values[...,power-1] += np.bincount(indices,
                                            weights=weight[part]**power, minlength=nx*np_).reshape(nx,np_)
                                    key = '__'.join((cat,tag,ch,REGIONS[ri]))
                                    entry = row['selected'].setdefault(key,dict(entries=0,sumw=0.,sumw2=0.,negative=0))
                                    entry['entries'] += int(part.sum())
                                    entry['sumw'] += float(weight[part].sum())
                                    entry['sumw2'] += float((weight[part]**2).sum())
                                    entry['negative'] += int((weight[part]<0).sum())
            if row['entries_read']!=tree.num_entries: raise ValueError('Incomplete input scan')
            audit.append(row)
            print(f"  done: {row['entries_read']} entries", flush=True)
    for cat in X: cache['mass_edges__'+cat] = X[cat]
    cache['pt_edges'] = PT
    np.savez_compressed(args.output/'source_histograms.npz', **cache)
    metadata = dict(date='2026-10-07', era='2018', channels=CHANNELS, regions=REGIONS,
                    processes=PROCESSES, pt_branches={'resolved':'pt2jet','merged':'ptmerged'},
                    config=config, config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest(),
                    producer_sha256=hashlib.sha256(Path(producer.__file__).read_bytes()).hexdigest(),
                    sf_file=str(sf_path), sf_sha256=hashlib.sha256(sf_path.read_bytes()).hexdigest(),
                    selections={'__'.join(key):v for key,v in expressions.items()},
                    samples=audit, array_axes=['channel','region','mass','pt','sumw_or_variance'],
                    ZH_subtracted=False, full_mass_tail=True)
    (args.output/'source_audit.json').write_text(json.dumps(metadata,indent=2))
    print('Source cache complete: '+str(args.output), flush=True)


if __name__=='__main__': main()
