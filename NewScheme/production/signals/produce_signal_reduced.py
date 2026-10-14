"""Signal reduced trees from separate DeepJet-augmented mass-specific ntuples."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[key]='1'
import argparse
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path
import time
import awkward as ak
import numpy as np
import uproot
import yaml
import background_helper_snapshot as helper

SAMPLES=[f'{mode}{mass}' for mode in ('ggh','vbf') for mass in (500,1000,3000)]
EXTRA=['weight','weight_no_pnet','weight_no_deepjet','weight_base','pnet_sf','pnet_sf_zvsqcd','pnet_sf_zvslight','deepjet_sf','w_deepjet','w_deepjet_up','w_deepjet_dn','w_deepjet_valid','nominal_mass','gen_mass','SumWeight','source_entry']

def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def validate(record,cat,tag,region,mass):
    assert all(np.isfinite(x).all() for x in record.values())
    np.testing.assert_array_equal(record['mass'],record['m2l2q'])
    assert (record['nominal_mass']==mass).all() and (record['SumWeight']>0).all()
    np.testing.assert_allclose(record['pnet_sf'],record['pnet_sf_zvsqcd']*record['pnet_sf_zvslight'],rtol=1e-14,atol=0)
    np.testing.assert_allclose(record['weight'],record['weight_no_pnet']*record['pnet_sf'],rtol=1e-13,atol=1e-14)
    np.testing.assert_allclose(record['weight_no_pnet'],record['weight_base']*record['deepjet_sf'],rtol=1e-13,atol=1e-14)
    np.testing.assert_allclose(record['weight_no_deepjet'],record['weight_base']*record['pnet_sf'],rtol=1e-13,atol=1e-14)
    if cat=='resolved' and tag!='vbftag':
        assert (record['w_deepjet_valid']==1).all() and (record['w_deepjet']>=0).all()
        np.testing.assert_array_equal(record['deepjet_sf'],record['w_deepjet'])
    else: assert (record['deepjet_sf']==1).all()
    if cat=='resolved': assert (record['pnet_sf']==1).all()
    if cat=='merged' and tag=='vbftag': assert (record['pnet_sf_zvslight']==1).all()
    q=record['mass2q']
    mask=((q>70)&(q<105)) if region=='SR' else (((q>135)&(q<165)) if region=='VR' else ((q<70)|(q>105)))
    assert mask.all(),'Hadronic region leakage'

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--era',choices=helper.ERAS,required=True); p.add_argument('--cat',choices=helper.CATEGORIES,required=True)
    p.add_argument('--input-dir',type=Path,required=True); p.add_argument('--output-dir',type=Path,required=True); p.add_argument('--sample',choices=SAMPLES)
    a=p.parse_args(); here=Path(__file__).resolve().parent
    cp=here/'cards'/f'config_UL{a.era}.yml'; config_hash=digest(cp); cfg=yaml.safe_load(cp.read_text())
    sfp=here/'sf'/f'NetSF_signal_{a.era}.yml'; sf_hash=digest(sfp); sfs=yaml.safe_load(sfp.read_text())
    samples=[a.sample] if a.sample else SAMPLES
    schema={key:'float64' for key in list(helper.BRANCHES[a.cat])+EXTRA}
    slots=[(sample,ch,tag,region) for sample in samples for ch in helper.CHANNELS for tag in helper.TAGS for region in helper.REGIONS]
    names={key:f'{key[0]}_{key[1]}_{a.cat}_{key[2]}_{a.era}_{key[3]}.root' for key in slots}
    a.output_dir.mkdir(parents=True,exist_ok=True)
    for name in names.values():
        if (a.output_dir/name).exists() or (a.output_dir/(name+'.partial')).exists(): raise FileExistsError(name)
    for ch in helper.CHANNELS.values():
        q='mass2jet' if a.cat=='resolved' else 'massmerged'
        for tag in helper.TAGS:
            assert f'(({q}<70) | ({q}>105))' in cfg['cut']['CR'][ch][a.cat][tag]
            assert f'{q}>135' in cfg['cut']['VR'][ch][a.cat][tag] and f'{q}<165' in cfg['cut']['VR'][ch][a.cat][tag]
    stats={k:dict(entries=0,sumw=0.,sumw2=0.,sumw_no_pnet=0.,bridge_entries=0,VR_entries=0,negative_weight_entries=0) for k in slots}
    input_records=[]; start=time.time()
    with ExitStack() as stack:
        handles={key:stack.enter_context(uproot.recreate(a.output_dir/(name+'.partial'))) for key,name in names.items()}
        trees={key:h.mktree('tree',schema) for key,h in handles.items()}
        for sample in samples:
            mass=int(sample[3:]); path=a.input_dir/a.era/f'{sample}.root'; audit=a.input_dir/a.era/f'{sample}.json'
            info=json.loads(audit.read_text()); assert info['status']=='passed' and info['nominal_mass']==mass
            xsec=float(cfg['samples_inf'][sample][1]); assert xsec>0
            fields=set(helper.BRANCHES[a.cat].values())|helper.selection_branches(cfg,helper.REGIONS,list(helper.CHANNELS),a.cat,helper.TAGS)
            fields.update(helper.CHANNELS.values()); fields.update(['EventWeight','SumWeight','GEN_H1_mass','w_deepjet','w_deepjet_up','w_deepjet_dn','w_deepjet_valid'])
            if a.cat=='merged': fields.update(['ptmerged','isbjet','iscjet','islightjet'])
            before=(path.stat().st_size,path.stat().st_mtime_ns)
            with uproot.open(path) as f:
                t=f['passedEvents']; assert not fields-set(t.keys())
                entry=0
                for events in t.iterate(sorted(fields),step_size=100000,library='ak'):
                    for ch,flag in helper.CHANNELS.items():
                        for region in helper.REGIONS:
                            membership=np.zeros(len(events),dtype=np.int8)
                            for tag in helper.TAGS:
                                expression=cfg['cut'][region][flag][a.cat][tag]
                                mask=ak.to_numpy(ak.numexpr.evaluate(expression,events)&events[flag]); membership+=mask
                                chosen=events[mask]
                                if not len(chosen): continue
                                base=ak.to_numpy(chosen['EventWeight']).astype(np.float64)*float(cfg['lumi'][a.era])*1000.*xsec/ak.to_numpy(chosen['SumWeight'])
                                dj=ak.to_numpy(chosen['w_deepjet']).astype(np.float64) if a.cat=='resolved' and tag!='vbftag' else np.ones(len(chosen))
                                qcd=helper.particle_net_sf(a.era,chosen,'ZvsQCD',sfs) if a.cat=='merged' else np.ones(len(chosen))
                                light=helper.particle_net_sf(a.era,chosen,'ZvsLight',sfs) if a.cat=='merged' and tag!='vbftag' else np.ones(len(chosen))
                                pnet=qcd*light
                                record={key:ak.to_numpy(chosen[field]).astype(np.float64) for key,field in helper.BRANCHES[a.cat].items()}
                                record.update(weight=base*dj*pnet,weight_no_pnet=base*dj,weight_no_deepjet=base*pnet,weight_base=base,
                                    pnet_sf=pnet,pnet_sf_zvsqcd=qcd,pnet_sf_zvslight=light,deepjet_sf=dj,nominal_mass=np.full(len(chosen),mass,dtype=np.float64),
                                    gen_mass=ak.to_numpy(chosen['GEN_H1_mass']).astype(np.float64),SumWeight=ak.to_numpy(chosen['SumWeight']).astype(np.float64),
                                    source_entry=(np.flatnonzero(mask)+entry).astype(np.float64))
                                for field in ('w_deepjet','w_deepjet_up','w_deepjet_dn','w_deepjet_valid'): record[field]=ak.to_numpy(chosen[field]).astype(np.float64)
                                validate(record,a.cat,tag,region,mass)
                                key=sample,ch,tag,region; trees[key].extend(record); row=stats[key]
                                row['entries']+=len(chosen); row['sumw']+=float(record['weight'].sum()); row['sumw2']+=float((record['weight']**2).sum()); row['sumw_no_pnet']+=float(record['weight_no_pnet'].sum())
                                q=record['mass2q']; row['bridge_entries']+=int(((q>105)&(q<135)).sum()); row['VR_entries']+=int(((q>135)&(q<165)).sum()); row['negative_weight_entries']+=int((record['weight']<0).sum())
                            assert (membership<=1).all(),'Tag overlap'
                    entry+=len(events)
            assert before==(path.stat().st_size,path.stat().st_mtime_ns)
            input_records.append(dict(sample=sample,path=str(path),entries=entry,size=before[0],mtime_ns=before[1],augmentation_report_sha256=digest(audit),xsec_pb=xsec))
            print(f'Completed {a.era} {a.cat} {sample}',flush=True)
    rows=[]
    for key,name in names.items():
        count=0; sums=np.zeros(3); previous=-1
        with uproot.open(a.output_dir/(name+'.partial')) as f:
            assert set(f['tree'].keys())==set(schema)
            for r in f['tree'].iterate(step_size=100000,library='np'):
                validate(r,a.cat,key[2],key[3],int(key[0][3:])); count+=len(r['weight'])
                assert (np.diff(r['source_entry'])>0).all()
                if len(r['source_entry']): assert r['source_entry'][0]>previous; previous=r['source_entry'][-1]
                sums+=np.array([r['weight'].sum(),(r['weight']**2).sum(),r['weight_no_pnet'].sum()])
        expected=stats[key]; assert count==expected['entries']
        np.testing.assert_allclose(sums,[expected[k] for k in ('sumw','sumw2','sumw_no_pnet')],rtol=1e-10,atol=1e-8)
        rows.append(dict(sample=key[0],channel=key[1],tag=key[2],region=key[3],file=name,**expected))
    assert digest(cp)==config_hash and digest(sfp)==sf_hash, 'Input configuration/calibration changed'
    report=dict(status='passed',era=a.era,category=a.cat,pilot=bool(a.sample),elapsed_seconds=time.time()-start,
        inputs=input_records,outputs=rows,branches=list(schema),config=str(cp),config_sha256=digest(cp),pnet_sf_sha256=digest(sfp),
        normalization='per-component per-event SumWeight; cross sections from mass-specific card entries (1 pb benchmark)',
        corrections='DeepJet nominal only resolved btag/untag; ParticleNet ZvsQCD all merged tags, ZvsLight merged btag/untag',
        broad_CR=True,VR_excluded_from_SR=True)
    label=f'{a.era}_{a.cat}'+(f'_{a.sample}_pilot' if a.sample else '')
    (a.output_dir/f'validation_{label}.json').write_text(json.dumps(report,indent=2))
    for name in names.values(): (a.output_dir/(name+'.partial')).rename(a.output_dir/name)
    (a.output_dir/f'{label}.done').touch()
    print(f'PASS {label}: {len(rows)} reduced trees fully read and validated',flush=True)
if __name__=='__main__': main()
