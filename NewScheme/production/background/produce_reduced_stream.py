"""Audited, bounded-memory all-background reduced trees; no alpha deployment."""
import argparse
from contextlib import ExitStack
import hashlib
import json
import logging
import os
from pathlib import Path
import time
os.environ.setdefault('NUMEXPR_NUM_THREADS','1')
import awkward as ak
import numpy as np
import uproot
import yaml
import make_reduced_trees as producer

SAMPLES = ['DY_pt50To100','DY_pt100To250','DY_pt250To400','DY_pt400To650',
           'DY_pt650ToInf','TTTo2L2Nu','WWTo2L2Nu','tZq','WZTo2Q2L','ZZTo2Q2L','Data']
WEIGHTS = ['weight','weight_no_pnet','pnet_sf','pnet_sf_zvsqcd','pnet_sf_zvslight']


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def required(config, category, sample):
    names = producer.selection_branches(config,producer.REGIONS,
        list(producer.CHANNELS),category,producer.TAGS)
    names.update(producer.BRANCHES[category].values())
    names.update(producer.CHANNELS.values())
    if sample != 'Data':
        names.add('EventWeight')
        if category=='resolved' and not sample.startswith('DY'):
            names.add('w_deepjet')
        if category=='merged' and producer.needs_particle_net(sample):
            names.update(['ptmerged','isbjet','iscjet','islightjet'])
    return sorted(names)


def preflight(config, args):
    records = []
    for sample in SAMPLES:
        path = producer.input_path(config,args.era,sample)
        if not path.is_file():
            raise FileNotFoundError(path)
        with uproot.open(path) as handle:
            tree = handle['passedEvents']
            names = required(config,args.cat,sample)
            missing = set(names)-set(tree.keys())
            if missing:
                raise KeyError(f'{sample}: missing {sorted(missing)}')
            norm = None if sample=='Data' else float(handle['sumWeights'].values().sum())
            if norm is not None and (not np.isfinite(norm) or norm==0):
                raise ValueError('Invalid sumWeights')
            records.append(dict(sample=sample,path=str(path),entries=tree.num_entries,
                sumWeights=norm,size=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns,branches=names))
    # Confirm the requested region convention, without changing other card cuts.
    for channel in producer.CHANNELS.values():
        mass = 'mass2jet' if args.cat=='resolved' else 'massmerged'
        for tag in producer.TAGS:
            cr = config['cut']['CR'][channel][args.cat][tag]
            vr = config['cut']['VR'][channel][args.cat][tag]
            assert f'(({mass}<70) | ({mass}>105))' in cr
            assert f'{mass}>135' in vr and f'{mass}<165' in vr
    return records


def slots(args):
    for process in ('DY','TTplusWW','VZ','Data'):
        for channel in producer.CHANNELS:
            for tag in producer.TAGS:
                for region in producer.REGIONS:
                    yield process,channel,tag,region


def validate(record, process, region):
    if any(not np.isfinite(a).all() for a in record.values()):
        raise ValueError('Nonfinite output branch')
    np.testing.assert_array_equal(record['mass'],record['m2l2q'])
    np.testing.assert_allclose(record['pnet_sf'],record['pnet_sf_zvsqcd']*record['pnet_sf_zvslight'],rtol=1e-14,atol=0)
    np.testing.assert_allclose(record['weight'],record['weight_no_pnet']*record['pnet_sf'],rtol=1e-14,atol=0)
    q = record['mass2q']
    if region=='CR': valid=(q<70)|(q>105)
    elif region=='SR': valid=(q>70)&(q<105)
    else: valid=(q>135)&(q<165)
    assert valid.all(), 'Hadronic region leakage'
    if process=='Data':
        for field in WEIGHTS:
            assert np.all(record[field]==1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--era',choices=producer.ERAS,required=True)
    parser.add_argument('--cat',choices=producer.CATEGORIES,required=True)
    parser.add_argument('--config-dir',type=Path,default=producer.DEFAULT_CONFIG_DIR)
    parser.add_argument('--sf-dir',type=Path,default=producer.DEFAULT_SF_DIR)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--step-size',type=int,default=100000)
    parser.add_argument('--max-entries',type=int)
    parser.add_argument('--preflight-only',action='store_true')
    args=parser.parse_args()
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(message)s')
    config_path=args.config_dir/f'config_UL{args.era}.yml'
    config=yaml.safe_load(config_path.read_text())
    inputs=preflight(config,args)
    sf_path=args.sf_dir/f'NetSF_signal_{args.era}.yml'
    sfs=yaml.safe_load(sf_path.read_text()) if args.cat=='merged' else None
    if args.preflight_only:
        print(json.dumps(dict(status='passed',era=args.era,category=args.cat,inputs=inputs),indent=2),flush=True)
        return
    out=args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    report_path=out/f'validation_{args.era}_{args.cat}.json'
    names={key:producer.output_name(key[0],key[1],args.cat,key[2],args.era,key[3]) for key in slots(args)}
    for name in names.values():
        if (out/name).exists() or (out/(name+'.partial')).exists():
            raise FileExistsError(f'Refusing to overwrite {name}')
    schema={field:'float64' for field in list(producer.BRANCHES[args.cat])+WEIGHTS}
    stats={key:dict(entries=0,sumw=0.,sumw2=0.,sumw_no_pnet=0.,bridge_entries=0,
                   negative_weight_entries=0) for key in names}
    start=time.time()
    with ExitStack() as stack:
        outputs={key:stack.enter_context(uproot.recreate(out/(name+'.partial'))) for key,name in names.items()}
        trees={key:handle.mktree('tree',schema) for key,handle in outputs.items()}
        for inp in inputs:
            sample=inp['sample']; group=producer.sample_group(sample)
            with uproot.open(inp['path']) as handle:
                stop=args.max_entries
                for events in handle['passedEvents'].iterate(inp['branches'],step_size=args.step_size,
                        entry_stop=stop,library='ak'):
                    for channel,ch_key in producer.CHANNELS.items():
                        for tag in producer.TAGS:
                            for region in producer.REGIONS:
                                expr=config['cut'][region][ch_key][args.cat][tag]
                                mask=ak.numexpr.evaluate(expr,events)&events[ch_key]
                                chosen=events[mask]
                                if not len(chosen): continue
                                weights=({key:np.ones(len(chosen)) for key in WEIGHTS} if sample=='Data' else
                                    producer.mc_weights(chosen,config,args.era,sample,group,args.cat,
                                        tag,inp['sumWeights'],sfs,return_components=True))
                                record={key:np.asarray(value,dtype=np.float64) for key,value in
                                    producer.reduced_record(chosen,args.cat,weights).items()}
                                validate(record,group,region)
                                key=group,channel,tag,region
                                trees[key].extend(record)
                                row=stats[key]; row['entries']+=len(chosen)
                                row['sumw']+=float(record['weight'].sum())
                                row['sumw2']+=float((record['weight']**2).sum())
                                row['sumw_no_pnet']+=float(record['weight_no_pnet'].sum())
                                row['bridge_entries']+=int(((record['mass2q']>105)&(record['mass2q']<135)).sum())
                                row['negative_weight_entries']+=int((record['weight']<0).sum())
            logging.info('Completed %s %s %s',args.era,args.cat,sample)
    rows=[]
    for key,name in names.items():
        expected=stats[key]; actual=np.zeros(3); n=0
        with uproot.open(out/(name+'.partial')) as handle:
            assert set(handle['tree'].keys())==set(schema)
            for record in handle['tree'].iterate(step_size=args.step_size,library='np'):
                validate(record,key[0],key[3]); n+=len(record['weight'])
                actual+=np.array([record['weight'].sum(),(record['weight']**2).sum(),record['weight_no_pnet'].sum()])
            assert n==expected['entries']
            np.testing.assert_allclose(actual,[expected['sumw'],expected['sumw2'],expected['sumw_no_pnet']],rtol=1e-11,atol=1e-7)
        rows.append(dict(file=name,process=key[0],channel=key[1],tag=key[2],region=key[3],**expected))
    for name in names.values():
        (out/(name+'.partial')).rename(out/name)
    report=dict(status='passed',era=args.era,category=args.cat,files=len(rows),schema=list(schema),
        sources=inputs,outputs=rows,elapsed_seconds=time.time()-start,pilot=args.max_entries is not None,
        config_sha256=digest(config_path),producer_sha256=digest(Path(producer.__file__)),
        runner_sha256=digest(__file__),particle_net_table=str(sf_path) if sfs is not None else None,
        particle_net_table_sha256=digest(sf_path) if sfs is not None else None,
        selections='Original era card cuts plus explicit channel guard; CR includes105-135 and VR; strict VR135<m2q<165. No m2l2q range cut.',
        calibration='Nominal era/flavour/pT/score SFs including DY; ZvsQCD for merged, additional ZvsLight for non-VBF merged.',
        weighting='weight=weight_no_pnet*pnet_sf; no-pnet preserves normalization, lepton and existing applicable DeepJet weights.',
        legacy_deepjet='Stored weights used unchanged, valid resolved non-VBF domain only. DY DeepJet remains unapplied.',
        alpha_outputs=False,SR_Data_plotted=False)
    report_path.write_text(json.dumps(report,indent=2))
    (out/f'{args.era}_{args.cat}.done').write_text('Validated '+str(len(rows))+' trees\n')
    print(json.dumps(dict(status='passed',era=args.era,category=args.cat,files=len(rows),
        entries=sum(row['entries'] for row in rows),elapsed_seconds=time.time()-start)),flush=True)


if __name__=='__main__':
    main()
