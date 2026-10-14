"""Produce only new ZHbb reduced trees, preserving all existing reductions."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
for name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[name] = '1'
import awkward as ak
import numpy as np
import uproot
import yaml


def region_mask(q, region):
    accepted = (q > 40) & (q < 180)
    sr = (q > 70) & (q < 105)
    vr = (q > 135) & (q < 165)
    if region == 'SR': return accepted & sr
    if region == 'VR': return accepted & vr
    return accepted & ((q < 70) | (q > 105)) & ~vr


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis', type=Path, default=Path('/cms/user/liugeliang/CMS/HZZAnalysis'))
    p.add_argument('--input', type=Path)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists(): raise FileExistsError('Use a new reduced-output directory')
    source = args.input or args.analysis/'ANATree/zhbb_2018/ZH_HToBB_ZToLL.root'
    sys.path.insert(0, str(args.analysis/'FitTemplate'))
    import make_reduced_trees as producer
    config_path = args.analysis/'cards/config_UL2018.yml'
    sf_path = producer.DEFAULT_SF_DIR/'NetSF_signal_2018.yml'
    config = yaml.safe_load(config_path.read_text())
    sfs = yaml.safe_load(sf_path.read_text())
    needed = {'EventWeight', 'SumWeight', 'w_deepjet', 'w_deepjet_valid',
              'isbjet', 'iscjet', 'islightjet', 'ptmerged'}
    for cat in producer.CATEGORIES:
        needed.update(producer.BRANCHES[cat].values())
        needed.update(producer.selection_branches(config, ['SR'], ['ee', 'mumu'], cat, producer.TAGS))
    with uproot.open(source) as f:
        missing = needed - set(f['passedEvents'].keys())
        if missing: raise ValueError('Missing branches: %s' % sorted(missing))
        events = f['passedEvents'].arrays(sorted(needed), library='ak')
        norm = float(f['sumWeights'].values().sum())
        processed = float(f['cutflow'].values()[0])
    if norm != 4547627 or processed != 4837847:
        raise ValueError('The expected validated full ZHbb sample was not supplied')
    if not np.all(ak.to_numpy(events['SumWeight']) == norm):
        raise ValueError('Per-event normalization disagrees with merged normalization')
    args.output.mkdir(parents=True)
    records = []
    for cat in producer.CATEGORIES:
        qname = producer.BRANCHES[cat]['mass2q']
        for channel in ['ee', 'mumu']:
            for tag in producer.TAGS:
                expression = config['cut']['SR'][producer.CHANNELS[channel]][cat][tag]
                expression, count = re.subn(
                    r'\(\('+qname+r'>70\)\s*&\s*\('+qname+r'<105\)\)',
                    '(('+qname+'>40) & ('+qname+'<180))', expression)
                if count != 1: raise ValueError('Unrecognized target mass selection')
                mask = ak.to_numpy(ak.numexpr.evaluate(expression, events))
                chosen = events[mask]
                unit = ak.to_numpy(chosen['EventWeight']).astype(float)*config['lumi']['2018']*1000*config['samples_inf']['ZH_HToBB'][1]/norm
                sf = np.ones(len(chosen))
                if cat == 'resolved' and tag != 'vbftag':
                    if not np.all(ak.to_numpy(chosen['w_deepjet_valid'])):
                        raise ValueError('Selected resolved event has invalid DeepJet weight')
                    sf *= ak.to_numpy(chosen['w_deepjet'])
                if cat == 'merged':
                    flavors = sum(ak.to_numpy(chosen[k]).astype(int) for k in ['isbjet', 'iscjet', 'islightjet'])
                    if not np.all(flavors == 1): raise ValueError('Merged flavour flags not one-hot')
                    sf *= producer.particle_net_sf('2018', chosen, 'ZvsQCD', sfs)
                    if tag != 'vbftag': sf *= producer.particle_net_sf('2018', chosen, 'ZvsLight', sfs)
                if not np.all(np.isfinite(sf)) or np.any(sf <= 0):
                    raise ValueError('Uncovered or invalid tag SF')
                q = ak.to_numpy(chosen[qname])
                for region in ['CR', 'SR', 'VR']:
                    selected = region_mask(q, region)
                    record = producer.reduced_record(chosen[selected], cat, (unit*sf)[selected])
                    record['weight_unit_tag'] = unit[selected]
                    path = args.output/producer.output_name('ZHbb', channel, cat, tag, '2018', region)
                    with uproot.recreate(path) as f: f['tree'] = record
                    records.append(dict(category=cat, channel=channel, tag=tag, region=region,
                        entries=int(selected.sum()), sumw=float(record['weight'].sum()),
                        unit_tag_sumw=float(record['weight_unit_tag'].sum()),
                        negative=int((record['weight'] < 0).sum()), path=str(path)))
                    print(records[-1], flush=True)
    metadata = dict(source=str(source), source_bytes=source.stat().st_size,
                    source_mtime=source.stat().st_mtime, normalization=norm, processed=processed,
                    lumi_fb=config['lumi']['2018'], xsec_pb=config['samples_inf']['ZH_HToBB'][1],
                    config=str(config_path), config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest(),
                    sf=str(sf_path), sf_sha256=hashlib.sha256(sf_path.read_bytes()).hexdigest(),
                    producer_sha256=hashlib.sha256(Path(producer.__file__).read_bytes()).hexdigest(),
                    cr='40<q<70 or 105<q<=135 or 165<=q<180; includes H window, excludes SR and VR',
                    merged_sf='Existing analysis ZvsQCD SF, plus ZvsLight for non-VBF; transfer to Hbb is a diagnostic convention, not an independently validated Hbb calibration',
                    records=records)
    (args.output/'metadata.json').write_text(json.dumps(metadata, indent=2))


if __name__ == '__main__': main()
