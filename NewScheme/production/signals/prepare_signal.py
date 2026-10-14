"""Non-destructive signal DeepJet augmentation with sample-specific efficiencies."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[key]='1'
import argparse
from array import array
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import correctionlib
import numpy as np
import uproot

BASE=Path(os.environ.get('HZZ_ANALYSIS_ROOT',Path(__file__).resolve().parents[3]))
SOURCE=Path('/cefs/higgs/zhuhao/HZZ_bkp/CMSSW_10_6_26/src/HZZAnalysis/UL_Legacy')
ERAS=('2016preAPV','2016postAPV','2017','2018')
WP=dict(zip(ERAS,(.2598,.2489,.3040,.2783)))
PT=np.array([30.,50.,70.,100.,140.,200.,300.,600.,1000.,np.inf])
ETA=np.array([0.,1.5,2.5])
FIELDS=['Met','foundZ1LCandidate','foundZ2JCandidate','isEE','isMuMu','KD_jjVBF','KD_jVBF']+[f'jet_{j}_{v}' for j in (1,2) for v in ('pt','eta','hadronflavor','deepbtag')]
VARIATIONS={'central':'w_deepjet','up':'w_deepjet_up','down':'w_deepjet_dn'}

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def domain(a):
    return (a['Met']<=150)&a['foundZ1LCandidate']&a['foundZ2JCandidate']&(a['isEE']|a['isMuMu'])&(a['KD_jjVBF']<.5)&(a['KD_jVBF']<.7)

def jet_sfs(a,mask,cal,variation):
    result=[]
    for j in (1,2):
        flavor=a[f'jet_{j}_hadronflavor'][mask].astype(np.int64)
        pt=a[f'jet_{j}_pt'][mask].astype(np.float64); eta=np.abs(a[f'jet_{j}_eta'][mask]).astype(np.float64)
        assert np.isin(flavor,[0,4,5]).all()
        out=np.empty(len(pt))
        for name,m in (('deepJet_comb',flavor!=0),('deepJet_incl',flavor==0)):
            if m.any(): out[m]=cal[name].evaluate(variation,'M',flavor[m],eta[m],pt[m])
        assert np.isfinite(out).all() and (out>0).all()
        result.append(out)
    return result

def check_reference(era,mode):
    path=SOURCE/era/'MC'/mode/'skimed'/f'{mode}.root'
    calpath=BASE/'deepjetSFsfiles'/era/'btagging.json.gz'
    cal=correctionlib.CorrectionSet.from_file(str(calpath))
    with uproot.open(path) as f:
        a=f['passedEvents'].arrays(FIELDS+list(VARIATIONS.values()),entry_stop=200000,library='np')
    selected=domain(a)&(a['jet_1_deepbtag']>=WP[era])&(a['jet_2_deepbtag']>=WP[era])
    assert selected.sum()>20,(era,mode,'insufficient reference b-tag events')
    differences={}
    for variation,name in VARIATIONS.items():
        sf=jet_sfs(a,selected,cal,variation); predicted=sf[0]*sf[1]; stored=a[name][selected]
        np.testing.assert_allclose(predicted,stored,rtol=2e-5,atol=1e-6)
        differences[name]=float(np.max(np.abs(predicted-stored)))
    return dict(era=era,mode=mode,path=str(path),btag_events=int(selected.sum()),max_abs_difference=differences,calibration_sha256=sha(calpath))

def clone_with_weights(source,target,weights):
    import ROOT
    if target.exists() or Path(str(target)+'.partial').exists(): raise FileExistsError(target)
    original=ROOT.TFile.Open(str(source),'READ'); tree=original.Get('passedEvents')
    assert tree and tree.GetEntries()==len(weights['w_deepjet'])
    for name in weights: assert not tree.GetBranch(name),name
    target.parent.mkdir(parents=True,exist_ok=True)
    out=ROOT.TFile(str(target)+'.partial','CREATE'); assert not out.IsZombie()
    out.cd(); clone=tree.CloneTree(-1,'fast')
    buffers={}; branches={}
    for name in weights:
        boolean=name=='w_deepjet_valid'; buffers[name]=array('b' if boolean else 'f',[0])
        branches[name]=clone.Branch(name,buffers[name],name+('/O' if boolean else '/F'))
    for i in range(tree.GetEntries()):
        for name,w in weights.items():
            buffers[name][0]=int(w[i]) if name=='w_deepjet_valid' else float(w[i]); branches[name].Fill()
    clone.Write('passedEvents',ROOT.TObject.kOverwrite)
    for key in original.GetListOfKeys():
        if key.GetName()!='passedEvents':
            obj=key.ReadObj(); out.cd(); obj.Write(key.GetName())
    out.Close(); original.Close()
    with uproot.open(source) as old,uproot.open(str(target)+'.partial') as new:
        assert new['passedEvents'].num_entries==old['passedEvents'].num_entries
        assert set(new['passedEvents'].keys())==set(old['passedEvents'].keys())|set(weights)
        for name,w in weights.items(): np.testing.assert_array_equal(new['passedEvents'][name].array(library='np'),w)
        for name in ('run','event','lumiSect','SumWeight','mass2l2jet','mass2lj','particleNetZvsQCD'):
            np.testing.assert_array_equal(new['passedEvents'][name].array(library='np'),old['passedEvents'][name].array(library='np'))
    Path(str(target)+'.partial').rename(target)

def prepare_one(task):
    era,mode,mass,output=task; out=Path(output)/era; sample=f'{mode}{mass}'
    paths=list((SOURCE/era/'MC'/mode/'skimed').glob(f'*_M{mass}_*.root')); assert len(paths)==1
    source=paths[0]; target=out/f'{sample}.root'
    if target.exists() or (out/f'{sample}.json').exists(): raise FileExistsError(target)
    before=(source.stat().st_size,source.stat().st_mtime_ns)
    with uproot.open(source) as f:
        a=f['passedEvents'].arrays(FIELDS,library='np'); n=f['passedEvents'].num_entries
        norm=f['passedEvents']['SumWeight'].array(library='np')
        assert np.isfinite(norm).all() and (norm>0).all() and (norm==norm[0]).all()
    selected=domain(a); assert selected.any()
    counts={}; passed_counts={}; efficiencies={}
    for flavor in (0,4,5):
        total=np.zeros((len(PT)-1,len(ETA)-1)); passed=total.copy()
        for j in (1,2):
            mask=selected&(a[f'jet_{j}_hadronflavor']==flavor)
            pt=a[f'jet_{j}_pt'][mask]; eta=np.abs(a[f'jet_{j}_eta'][mask]); tagged=a[f'jet_{j}_deepbtag'][mask]>=WP[era]
            total+=np.histogram2d(pt,eta,[PT,ETA])[0]; passed+=np.histogram2d(pt[tagged],eta[tagged],[PT,ETA])[0]
        counts[flavor]=total; passed_counts[flavor]=passed
        efficiencies[flavor]=np.divide(passed,total,out=np.full_like(total,np.nan),where=total>0)
    eps=[]; passes=[]
    for j in (1,2):
        pt=a[f'jet_{j}_pt'][selected]; eta=np.abs(a[f'jet_{j}_eta'][selected]); flavor=a[f'jet_{j}_hadronflavor'][selected]
        ix=np.searchsorted(PT,pt,side='right')-1; iy=np.searchsorted(ETA,eta,side='right')-1
        assert ((ix>=0)&(ix<len(PT)-1)&(iy>=0)&(iy<len(ETA)-1)).all()
        e=np.empty(len(pt))
        for f in (0,4,5):
            mask=flavor==f; e[mask]=efficiencies[f][ix[mask],iy[mask]]
        assert np.isfinite(e).all()
        passed=a[f'jet_{j}_deepbtag'][selected]>=WP[era]
        assert (e[~passed]<1).all()
        eps.append(e); passes.append(passed)
    calpath=BASE/'deepjetSFsfiles'/era/'btagging.json.gz'; cal=correctionlib.CorrectionSet.from_file(str(calpath))
    weights={}; stats={}
    for variation,name in VARIATIONS.items():
        sfs=jet_sfs(a,selected,cal,variation); jetweights=[]
        for sf,e,passed in zip(sfs,eps,passes):
            w=sf.copy(); failed=~passed; w[failed]=(1-sf[failed]*e[failed])/(1-e[failed]); jetweights.append(w)
        w=jetweights[0]*jetweights[1]; assert np.isfinite(w).all()
        if variation=='central': assert (w>=0).all(),(sample,era,'negative central probability weight')
        full=np.full(n,-999.,dtype=np.float32); full[selected]=w.astype(np.float32); weights[name]=full
        stats[name]=dict(min=float(w.min()),max=float(w.max()),mean=float(w.mean()),negative=int((w<0).sum()))
    weights['w_deepjet_valid']=selected
    out.mkdir(parents=True,exist_ok=True)
    efficiency={'pt_edges':PT,'eta_edges':ETA}
    for f in counts:
        for name,data in (('total',counts),('pass',passed_counts),('epsilon',efficiencies)): efficiency[f'{name}_{f}']=data[f]
    np.savez_compressed(out/f'{sample}_efficiency.npz',**efficiency)
    clone_with_weights(source,target,weights)
    assert before==(source.stat().st_size,source.stat().st_mtime_ns),'Original ntuple changed'
    report=dict(status='passed',era=era,sample=sample,nominal_mass=mass,source=str(source),output=str(target),entries=n,
        selected=int(selected.sum()),SumWeight=float(norm[0]),WP=WP[era],calibration=str(calpath),calibration_sha256=sha(calpath),
        weights=stats,efficiency='unweighted pass/total jets in this mass-specific signal sample, nominal resolved/non-VBF domain; no borrowed efficiency or clipping',
        source_size=before[0],source_mtime_ns=before[1],event_order_and_unrelated_branches_preserved=True)
    (out/f'{sample}.json').write_text(json.dumps(report,indent=2))
    print(f'PASS augmented {era} {sample}: {n} entries, {int(selected.sum())} DeepJet-valid',flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--output-dir',type=Path,required=True); p.add_argument('--pilot',action='store_true'); p.add_argument('--reference-only',action='store_true')
    args=p.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True)
    if args.reference_only:
        reports=[check_reference(era,mode) for era in ERAS for mode in ('ggh','vbf')]
        (args.output_dir/'deepjet_reference_validation.json').write_text(json.dumps(reports,indent=2)); print('PASS: all-era DeepJet calibration regression',flush=True); return
    tasks=[(era,mode,mass,str(args.output_dir)) for era in ERAS for mode in ('ggh','vbf') for mass in (500,1000,3000)]
    if args.pilot: tasks=[('2018','ggh',500,str(args.output_dir))]
    with ProcessPoolExecutor(max_workers=2) as executor: reports=list(executor.map(prepare_one,tasks))
    (args.output_dir/'augmentation_manifest.json').write_text(json.dumps(reports,indent=2))
    print(f'Completed {len(reports)} non-destructive augmented signal ntuples',flush=True)
if __name__=='__main__': main()
