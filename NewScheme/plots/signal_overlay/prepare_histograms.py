"""Cached Run-2 signal overlays and blinded direct-MC mass spectra."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import argparse,csv,hashlib,json,re
from pathlib import Path
import numpy as np
import numexpr as ne
import uproot,yaml

ERAS=('2016preAPV','2016postAPV','2017','2018')
TAGS=('untag','btag','vbftag'); CHANNELS=('ee','mumu'); MC=('VZ','TTplusWW','DY')
SIGNALS=tuple(f'{p}{m}' for p in ('ggh','vbf') for m in (500,1000,3000))
VARS=('kd','kd_jjvbf','kd_jvbf','particleNetZvsQCD','particleNetZbbvslight')
SCORE_EDGES=np.linspace(0,1,21)
MASS_EDGES=np.r_[np.arange(300.,1500.,50.),np.arange(1500.,3501.,100.),np.inf]

def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def fill(store,key,x,w,edges):
    if key not in store: store[key]=np.zeros((2,len(edges)-1))
    store[key][0]+=np.histogram(x,edges,weights=w)[0]
    store[key][1]+=np.histogram(x,edges,weights=w*w)[0]

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--analysis',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
    sigroot=a.analysis/'FitTemplate'/'RedTrees_signal_20261012'; bgroot=a.analysis/'FitTemplate'/'RedTrees_20261009_pnet'
    archive={}; mass={}; full={}; audits=[]; rows=[]; hashes={}
    for root,label in ((sigroot,'signal'),(bgroot,'background')):
        for era in ERAS:
            for cat in ('resolved','merged'):
                path=root/f'validation_{era}_{cat}.json'; report=json.loads(path.read_text())
                assert report['status']=='passed' and not report['pilot'] and (root/f'{era}_{cat}.done').is_file()
                hashes[f'{label}_{path.name}']=digest(path)
    for cat in ('resolved','merged'):
        for region in ('SR','CR'):
            for era in ERAS:
                for tag in TAGS:
                    for ch in CHANNELS:
                        for sample in SIGNALS:
                            path=sigroot/f'{sample}_{ch}_{cat}_{tag}_{era}_{region}.root'
                            audits.append(dict(kind='signal_reduced',file=str(path),size=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
                            n=0; sums=np.zeros(2)
                            with uproot.open(path) as f:
                                for v in f['tree'].iterate(list(VARS)+['mass2q','m2l2q','weight','weight_no_pnet'],step_size=100000,library='np'):
                                    q=v['mass2q']; selected=(v['m2l2q']>=300)&(q>=40)&(q<180)
                                    w=v['weight'][selected].astype(np.float64); no=v['weight_no_pnet'][selected].astype(np.float64)
                                    n+=int(selected.sum()); sums+=np.array([w.sum(),(w*w).sum()])
                                    assert np.isfinite(w).all()
                                    for var in VARS:
                                        x=v[var][selected]; good=np.isfinite(x)&(x>=0)&(x<=1)
                                        fill(archive,f'{cat}_{region}_{var}__{sample}',x[good],w[good],SCORE_EDGES)
                                        if cat=='merged' and var.startswith('particleNet'):
                                            fill(archive,f'{cat}_{region}_{var}__{sample}_no_pnet',x[good],no[good],SCORE_EDGES)
                                    if region=='SR': fill(mass,f'{cat}__{tag}__{ch}__{sample}',v['m2l2q'][selected],w,MASS_EDGES)
                            rows.append(dict(kind='signal_reduced',era=era,category=cat,region=region,tag=tag,channel=ch,sample=sample,entries=n,sumw=sums[0],sumw2=sums[1]))
            print('Loaded reduced signals',cat,region,flush=True)
    for cat in ('resolved','merged'):
        for era in ERAS:
            for tag in TAGS:
                for ch in CHANNELS:
                    for sample in MC:
                        path=bgroot/f'{sample}_{ch}_{cat}_{tag}_{era}_SR.root'
                        audits.append(dict(kind='background_SR_MC',file=str(path),size=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
                        n=0; sums=np.zeros(2)
                        with uproot.open(path) as f:
                            for v in f['tree'].iterate(['mass2q','m2l2q','weight'],step_size=100000,library='np'):
                                q=v['mass2q']; assert ((q>70)&(q<105)).all()
                                selected=(v['m2l2q']>=300)&(q>=40)&(q<180); w=v['weight'][selected].astype(np.float64)
                                assert np.isfinite(w).all()
                                fill(mass,f'{cat}__{tag}__{ch}__{sample}',v['m2l2q'][selected],w,MASS_EDGES)
                                n+=int(selected.sum()); sums+=np.array([w.sum(),(w*w).sum()])
                        rows.append(dict(kind='background_SR_MC',era=era,category=cat,region='SR',tag=tag,channel=ch,sample=sample,entries=n,sumw=sums[0],sumw2=sums[1]))
        print('Loaded SR background masses',cat,flush=True)
    for era in ERAS:
        cp=a.analysis/'NewScheme'/'production'/'signals'/'cards'/f'config_UL{era}.yml'; cfg=yaml.safe_load(cp.read_text()); hashes[cp.name]=digest(cp)
        cuts={}
        for region in ('SR','CR'):
            for ch,flag in (('ee','isEE'),('mumu','isMuMu')):
                for tag in TAGS:
                    expr=cfg['cut'][region][flag]['merged'][tag]
                    expr,n=re.subn(r'\(particleNetZvsQCD\s*>\s*0\.9\)','(True)',expr); assert n==1
                    cuts[(region,ch,tag)]=f'({expr}) & ({flag}==True)'
        fields={'particleNetZvsQCD','massmerged','mass2lj','EventWeight','SumWeight'}
        for expr in cuts.values(): fields.update(re.findall(r'[A-Za-z_][A-Za-z0-9_]*',expr))
        fields.difference_update({'True','False'})
        for sample in SIGNALS:
            path=a.analysis/'FitTemplate'/'SignalWithDeepJet_20261012'/era/f'{sample}.root'
            audits.append(dict(kind='signal_fullrange',file=str(path),size=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
            counts={region:np.zeros(3) for region in ('SR','CR')}
            with uproot.open(path) as f:
                for v in f['passedEvents'].iterate(sorted(fields),step_size=100000,library='np'):
                    x=v['particleNetZvsQCD']; q=v['massmerged']
                    baseline=(v['mass2lj']>=300)&(q>=40)&(q<180)&np.isfinite(x)&(x>=0)&(x<=1)
                    w=v['EventWeight'].astype(np.float64)*float(cfg['lumi'][era])*1000.*float(cfg['samples_inf'][sample][1])/v['SumWeight']
                    assert np.isfinite(w[baseline]).all()
                    for region in ('SR','CR'):
                        member=np.zeros(len(x),dtype=np.int8)
                        for key,expr in cuts.items():
                            if key[0]==region: member+=ne.evaluate(expr,local_dict=v).astype(np.int8)
                        assert (member<=1).all()
                        mask=baseline&(member==1); prefix=f'merged_{region}_particleNetZvsQCD__{sample}'
                        fill(full,prefix,x[mask],w[mask],SCORE_EDGES)
                        passed=mask&(x>.9); fill(full,prefix+'__tagged_check',x[passed],w[passed],SCORE_EDGES)
                        counts[region]+=np.array([mask.sum(),w[mask].sum(),(w[mask]**2).sum()])
            for region,c in counts.items(): rows.append(dict(kind='signal_fullrange',era=era,category='merged',region=region,tag='all',channel='combined',sample=sample,entries=int(c[0]),sumw=c[1],sumw2=c[2]))
        print('Loaded full-range signals',era,flush=True)
    for region in ('SR','CR'):
        for sample in SIGNALS:
            prefix=f'merged_{region}_particleNetZvsQCD__{sample}'
            np.testing.assert_allclose(full[prefix+'__tagged_check'],archive[prefix+'_no_pnet'],rtol=1e-10,atol=1e-7)
    for cat in ('resolved','merged'):
        for tag in TAGS:
            for ch in CHANNELS:
                for sample in MC+SIGNALS:
                    wanted=[r for r in rows if (r['category'],r['region'],r['tag'],r['channel'],r['sample'])==(cat,'SR',tag,ch,sample) and r['kind']!='signal_fullrange']
                    np.testing.assert_allclose(mass[f'{cat}__{tag}__{ch}__{sample}'].sum(axis=1),[sum(r[k] for r in wanted) for k in ('sumw','sumw2')],rtol=1e-10,atol=1e-7)
    assert not any('Data_' in r['file'] and '_SR.root' in r['file'] for r in audits)
    np.savez_compressed(a.output/'signal_histograms.npz',**archive)
    np.savez_compressed(a.output/'signal_fullrange_histograms.npz',**full)
    np.savez_compressed(a.output/'mass_histograms.npz',mass_edges=MASS_EDGES,**mass)
    with (a.output/'input_yields.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    (a.output/'metadata.json').write_text(json.dumps(dict(source_files=audits,source_hashes=hashes,SR_Data_read=False,signals=list(SIGNALS),signal_xsec_pb=1,
        style='ggF dotted; VBF solid; 500 red, 1000 dark red, 3000 pink; same mass colors for both modes',
        normalization='original era luminosities, displayed as138 fb^-1; full signal nominal weights; full-range no ParticleNet SF',
        mass_binning='common event-level bins:50 GeV300-1500,100 GeV1500-3500,>=3500 overflow displayed3500-4000; Events/50GeV',
        mass_background='direct DY MC with nominal ParticleNet SF, TTplusWW and VZ; no alpha prediction'),indent=2))
    print('PASS: SR blinding, full-range/tagged signal regression, full mass-spectrum yields and variances',flush=True)
if __name__=='__main__': main()
