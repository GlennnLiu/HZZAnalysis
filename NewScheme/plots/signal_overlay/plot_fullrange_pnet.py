"""Full-range pre-working-point ParticleNet ZvsQCD, original ntuples, no PNet SF."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
import argparse
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import mplhep as hep
import numpy as np
from scipy.stats import chi2
import uproot
import yaml
import numexpr as ne
import re
from concurrent.futures import ProcessPoolExecutor

ERAS=('2016preAPV','2016postAPV','2017','2018')
TAGS=('untag','btag','vbftag')
MC=('VZ','TTplusWW','DY')
VARIABLES=('particleNetZvsQCD',)
COLORS={'VZ':'#f4b183','TTplusWW':'#b1a2c5','DY':'#80c17d'}
LABELS={'VZ':'VZ','TTplusWW':r'$t\bar{t}$ + WW + tZq','DY':'DY + jets'}
AXES={'particleNetZvsQCD':'ParticleNet Z vs QCD'}
EDGES=np.linspace(0.,1.,21)

def group_of(sample):
    if sample=='Data': return 'Data'
    if sample.startswith('DY'): return 'DY'
    if sample in ('TTTo2L2Nu','WWTo2L2Nu','tZq'): return 'TTplusWW'
    if sample in ('WZTo2Q2L','ZZTo2Q2L'): return 'VZ'
    raise ValueError(sample)

def read_era(task):
    era,config_path=task
    text=Path(config_path).read_text(); cfg=yaml.safe_load(text)
    score='particleNetZvsQCD'
    cuts={}
    for region in ('SR','CR'):
        for channel,flag in (('ee','isEE'),('mumu','isMuMu')):
            for tag in TAGS:
                old=cfg['cut'][region][flag]['merged'][tag]
                new,n=re.subn(r'\(particleNetZvsQCD\s*>\s*0\.9\)','(True)',old)
                assert n==1, (era,region,channel,tag,old)
                cuts[(region,channel,tag)]=f'({new}) & ({flag}==True)'
    archive={}
    for region in ('SR','CR'):
        for proc in MC+(() if region=='SR' else ('Data',)):
            key=f'merged_{region}_{score}__{proc}'
            archive[key]=np.zeros((2,20))
            archive[key+'__tagged_check']=np.zeros((2,20))
    rows=[]; files=[]
    # Match the validated streaming production's explicit sample inventory;
    # do not change Samples_lists in any original YAML.
    samples=('DY_pt50To100','DY_pt100To250','DY_pt250To400','DY_pt400To650',
             'DY_pt650ToInf','TTTo2L2Nu','WWTo2L2Nu','tZq','WZTo2Q2L','ZZTo2Q2L','Data')
    for sample in samples:
        group=group_of(sample)
        regions=('CR',) if sample=='Data' else ('SR','CR')
        expressions={key:value for key,value in cuts.items() if key[0] in regions}
        fields={'massmerged','mass2lj',score}
        for expression in expressions.values():
            fields.update(re.findall(r'[A-Za-z_][A-Za-z0-9_]*',expression))
        fields.difference_update({'True','False'})
        if sample!='Data': fields.add('EventWeight')
        root=Path(cfg['ori_dir'])/era
        path=(root/'Data'/'skimed'/f'Data{era}UL_noDuplicates.root') if sample=='Data' else (root/cfg['samples_inf'][sample][0]/'skimed'/f'{sample}.root')
        files.append(dict(era=era,sample=sample,file=str(path),size=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
        stats={region:dict(selected=0,valid=0,invalid=0,below_wp=0,tagged=0,sumw=0.,sumw2=0.,tagged_sumw=0.,tagged_sumw2=0.,bridge_entries=0,VR_entries=0) for region in regions}
        with uproot.open(path) as f:
            normalization=1. if sample=='Data' else float(cfg['lumi'][era])*1000.*float(cfg['samples_inf'][sample][1])/float(f['sumWeights'].to_boost().sum())
            assert np.isfinite(normalization) and normalization>0
            tree=f['passedEvents']
            missing=fields-set(tree.keys()); assert not missing,(path,missing)
            for a in tree.iterate(sorted(fields),step_size=200000,library='np'):
                x=a[score]; q=a['massmerged']
                baseline=(a['mass2lj']>=300)&(q>=40)&(q<180)
                weights=np.ones(len(x),dtype=np.float64) if sample=='Data' else a['EventWeight'].astype(np.float64)*normalization
                for region in regions:
                    membership=np.zeros(len(x),dtype=np.int8)
                    for key,expression in expressions.items():
                        if key[0]==region: membership+=ne.evaluate(expression,local_dict=a).astype(np.int8)
                    assert (membership<=1).all(),(era,sample,region,'overlapping tag/channel cuts')
                    selected=baseline&(membership==1)
                    c=stats[region]; c['selected']+=int(selected.sum())
                    valid=selected&np.isfinite(x)&(x>=0)&(x<=1)
                    assert np.isfinite(weights[valid]).all()
                    c['valid']+=int(valid.sum()); c['invalid']+=int((selected&~valid).sum())
                    passing=valid&(x>.9); c['tagged']+=int(passing.sum()); c['below_wp']+=int((valid&(x<=.9)).sum())
                    for mask,suffix in ((valid,''),(passing,'__tagged_check')):
                        w=weights[mask]; key=f'merged_{region}_{score}__{group}'+suffix
                        archive[key][0]+=np.histogram(x[mask],EDGES,weights=w)[0]
                        archive[key][1]+=np.histogram(x[mask],EDGES,weights=w*w)[0]
                        pre='tagged_' if suffix else ''
                        c[pre+'sumw']+=float(w.sum()); c[pre+'sumw2']+=float((w*w).sum())
                    c['bridge_entries']+=int((valid&(q>105)&(q<135)).sum())
                    c['VR_entries']+=int((valid&(q>135)&(q<165)).sum())
        for region,c in stats.items():
            assert c['selected']==c['valid']+c['invalid'] and c['valid']==c['tagged']+c['below_wp']
            rows.append(dict(era=era,sample=sample,process=group,region=region,normalization=normalization,**c))
        print(f'Processed {era} {sample}',flush=True)
    config=dict(era=era,path=str(config_path),sha256=hashlib.sha256(text.encode()).hexdigest(),original_yaml=text,cuts=[dict(region=k[0],channel=k[1],tag=k[2],expression=v) for k,v in cuts.items()])
    # Retain completed per-era results before final regression checks.
    cache=Path(__file__).resolve().parents[3]/'Plot'/'run2_PNet_fullrange'/f'cache_{era}'
    np.savez_compressed(str(cache)+'.npz',**archive)
    Path(str(cache)+'.json').write_text(json.dumps(dict(rows=rows,files=files,config=config)))
    return archive,rows,files,config

def load(source,out):
    tasks=[(era,str(source/'FitTemplate'/'jobs'/'reduced_20261009_pnet'/'cards'/f'config_UL{era}.yml')) for era in ERAS]
    archive={}; rows=[]; files=[]; configs=[]
    with ProcessPoolExecutor(max_workers=2) as executor:
        for result,era_rows,era_files,config in executor.map(read_era,tasks):
            for key,h in result.items(): archive[key]=archive.get(key,np.zeros_like(h))+h
            rows+=era_rows; files+=era_files; configs.append(config)
    for region in ('SR','CR'):
        prefix=f'merged_{region}_particleNetZvsQCD'
        archive[prefix+'__edges']=EDGES.copy()
        for proc in MC+(() if region=='SR' else ('Data',)):
            subset=[r for r in rows if r['region']==region and r['process']==proc]
            for suffix,pre in (('',''),('__tagged_check','tagged_')):
                np.testing.assert_allclose(archive[prefix+'__'+proc+suffix].sum(axis=1),[sum(r[pre+k] for r in subset) for k in ('sumw','sumw2')],rtol=1e-10,atol=1e-6)
    reference=source/'Plot'/'run2_PNet_scores'/'histograms.npz'
    with np.load(reference) as prior:
        for region in ('SR','CR'):
            prefix=f'merged_{region}_particleNetZvsQCD'
            np.testing.assert_allclose(archive[prefix+'__DY__tagged_check'],prior[prefix+'__DY_no_pnet'],rtol=3e-6,atol=1e-5)
        np.testing.assert_array_equal(archive['merged_CR_particleNetZvsQCD__Data__tagged_check'],prior['merged_CR_particleNetZvsQCD__Data'])
    assert not any('_SR_' in k and '__Data' in k for k in archive)
    np.savez_compressed(out/'histograms.npz',**archive)
    with (out/'input_yields.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    (out/'metadata.json').write_text(json.dumps(dict(source_files=files,config_snapshots=configs,luminosity_display=138,SR_Data_used=False,
        weighting='EventWeight*lumi*1000*xsec/sumWeights; no ParticleNet SF for any MC',
        selection='original merged tag/channel cuts with only particleNetZvsQCD>0.9 removed; m2l2q>=300 and 40<=massmerged<180; full CR includes bridge and VR',
        tagged_validation='raw DY passing region matches prior no-PNet DY bin contents/variances; CR Data passing region matches prior Data exactly',
        note='Original Data ntuples contain all regions; only CR selections are evaluated and histogrammed; no SR Data counts retained'),indent=2))
    print('PASS: yields/variances, SR blinding, and tagged-region agreement with prior plots',flush=True)
    return archive

def poisson(n):
    lo=np.zeros_like(n); mask=n>0; lo[mask]=.5*chi2.ppf(.15865525393145707,2*n[mask])
    hi=.5*chi2.ppf(.8413447460685429,2*(n+1))
    return np.stack((n-lo,hi-n))

def draw(archive,cat,region,var,axes=None):
    if axes is None:
        if region=='CR':
            fig,(ax,ratio)=plt.subplots(2,1,figsize=(11,10),sharex=True,gridspec_kw={'height_ratios':[3.3,1],'hspace':.05})
        else: fig,ax=plt.subplots(figsize=(11,9)); ratio=None
        fig.subplots_adjust(left=.18,right=.97,top=.90,bottom=.16)
    else: ax,ratio=axes; fig=ax.figure
    total=np.zeros(20); variance=np.zeros(20)
    for proc in MC:
        h=archive[f'{cat}_{region}_{var}__{proc}']
        ax.bar(EDGES[:-1],h[0],width=np.diff(EDGES),bottom=total,align='edge',color=COLORS[proc],edgecolor='none',label=LABELS[proc])
        total+=h[0]; variance+=h[1]
    err=np.sqrt(variance); step=lambda a:np.r_[a,a[-1]]
    band=r'Bkg. MC stat. ($1\sigma$)'
    ax.fill_between(EDGES,step(total-err),step(total+err),step='post',color='#4c566a',alpha=.3,lw=0,label=band)
    maximum=float(np.max(total+err)); centers=(EDGES[:-1]+EDGES[1:])/2
    if region=='CR':
        data=archive[f'{cat}_{region}_{var}__Data'][0]; de=poisson(data)
        ax.errorbar(centers,data,yerr=de,fmt='o',color='black',ms=5,lw=1.2,label='Data',zorder=5)
        maximum=max(maximum,float((data+de[1]).max())); good=total>0
        rel=np.divide(err,total,out=np.full_like(err,np.nan),where=good)
        ratio.fill_between(EDGES,step(1-rel),step(1+rel),step='post',color='#4c566a',alpha=.3,lw=0)
        ratio.errorbar(centers[good],data[good]/total[good],yerr=de[:,good]/total[good],fmt='o',color='black',ms=5,lw=1.2)
        off=good&(data>2.5*total); ratio.plot(centers[off],np.full(int(off.sum()),2.43),'^',color='black',ms=7)
        ratio.axhline(1,color='.4',lw=1); ratio.set_ylim(0,2.5); ratio.set_ylabel('Data / MC',fontsize=27)
    positive=total[total>0]
    ax.set_yscale('log'); ax.set_ylim(0.1,max(maximum,1.)*10000)
    hep.cms.label('Preliminary',data=region=='CR',lumi=138,com=13,ax=ax,fontsize=27)
    ax.text(.035,.95,f'Run 2, {cat.capitalize()}\nAll tags, '+r'$ee+\mu\mu$'+'\n'+('SR (blinded)' if region=='SR' else 'CR')+'\nNo PNet SF',transform=ax.transAxes,va='top',fontsize=23)
    handles,labels=ax.get_legend_handles_labels(); lookup=dict(zip(labels,handles)); order=[LABELS[p] for p in MC]+[band]+(['Data'] if region=='CR' else [])
    ax.legend([lookup[k] for k in order],order,loc='upper right',fontsize=22,frameon=False,handlelength=1.2,labelspacing=.3)
    ax.set_ylabel('Events / 0.05',fontsize=32)
    bottom=ratio if ratio is not None else ax
    bottom.set_xlabel(AXES[var],fontsize=34); bottom.set_xlim(0,1); bottom.set_xticks(np.linspace(0,1,6))
    for axis in (ax,ratio):
        if axis is not None: axis.tick_params(direction='in',top=True,right=True,labelsize=23); axis.grid(False)
    if ratio is not None: ax.tick_params(labelbottom=False)
    return fig

def render(archive,out):
    plt.style.use(hep.style.CMS)
    var='particleNetZvsQCD'; cat='merged'
    with PdfPages(out/'all_Run2_fullrange_PNet.pdf') as bundle:
        fig=plt.figure(figsize=(22,10)); outer=fig.add_gridspec(1,2,left=.08,right=.985,bottom=.16,top=.90,wspace=.32)
        for i,region in enumerate(('SR','CR')):
            one=draw(archive,cat,region,var)
            for ext in ('png','pdf'): one.savefig(out/f'merged_{region}_particleNetZvsQCD_fullrange.{ext}',dpi=140)
            bundle.savefig(one); plt.close(one)
            if region=='CR':
                sub=outer[i].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05); ax=fig.add_subplot(sub[0]); ratio=fig.add_subplot(sub[1],sharex=ax)
            else: ax=fig.add_subplot(outer[i]); ratio=None
            draw(archive,cat,region,var,(ax,ratio))
        for ext in ('png','pdf'): fig.savefig(out/f'SR_CR_fullrange_particleNetZvsQCD.{ext}',dpi=120)
        plt.close(fig)
    print('Completed full-range SR/CR plots and overview',flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--input',type=Path,required=True); p.add_argument('--output',type=Path,required=True); p.add_argument('--histograms',type=Path)
    a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
    if a.histograms:
        with np.load(a.histograms) as old: archive={k:old[k].copy() for k in old.files}
    else: archive=load(a.input,a.output)
    render(archive,a.output)
if __name__=='__main__': main()
