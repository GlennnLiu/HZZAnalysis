"""Run-2 ParticleNet scores, DY SF on/off, SR Data never opened."""
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

ERAS=('2016preAPV','2016postAPV','2017','2018')
TAGS=('untag','btag','vbftag')
MC=('VZ','TTplusWW','DY')
VARIABLES=('particleNetZvsQCD','particleNetZbbvslight')
COLORS={'VZ':'#f4b183','TTplusWW':'#b1a2c5','DY':'#80c17d'}
LABELS={'VZ':'VZ','TTplusWW':r'$t\bar{t}$ + WW + tZq','DY':'DY + jets, SF on'}
AXES={'particleNetZvsQCD':'ParticleNet Z vs QCD', 'particleNetZbbvslight':r'ParticleNet $Z\to b\bar{b}$ vs light'}
EDGES=np.linspace(0.,1.,21)

def load(source,out):
    archive={}; rows=[]; sources=[]; hashes={}
    for cat in ('merged',):
        for era in ERAS:
            p=source/f'validation_{era}_{cat}.json'; report=json.loads(p.read_text())
            assert report['status']=='passed' and not report['pilot'] and (source/f'{era}_{cat}.done').is_file()
            hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
        for region in ('SR','CR'):
            processes=MC if region=='SR' else MC+('Data',)
            for var in VARIABLES:
                prefix=f'{cat}_{region}_{var}'; archive[prefix+'__edges']=EDGES.copy()
                for proc in processes+('DY_no_pnet',): archive[prefix+'__'+proc]=np.zeros((2,20))
            for era in ERAS:
                for tag in TAGS:
                    for channel in ('ee','mumu'):
                        for proc in processes:
                            path=source/f'{proc}_{channel}_{cat}_{tag}_{era}_{region}.root'
                            sources.append(dict(file=str(path),size=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
                            counts={v:dict(selected=0,valid=0,negative=0,nonfinite=0,above_one=0,sumw=0.,sumw2=0.) for v in VARIABLES}
                            off_counts={v:dict(sumw=0.,sumw2=0.) for v in VARIABLES}
                            with uproot.open(path) as f:
                                for a in f['tree'].iterate(list(VARIABLES)+['mass2q','m2l2q','weight','weight_no_pnet','pnet_sf'],step_size=100000,library='np'):
                                    q=a['mass2q']; select=(a['m2l2q']>=300)&(q>=40)&(q<180)
                                    assert (((q>70)&(q<105)) if region=='SR' else ((q<70)|(q>105))).all()
                                    w=a['weight'][select]; assert np.isfinite(w).all()
                                    if proc=='Data': np.testing.assert_array_equal(w,np.ones(len(w)))
                                    wo=a['weight_no_pnet'][select] if proc=='DY' else w
                                    assert np.isfinite(wo).all()
                                    if proc=='DY': np.testing.assert_allclose(w,wo*a['pnet_sf'][select],rtol=2e-6,atol=1e-7)
                                    good={}
                                    for var in VARIABLES:
                                        x=a[var][select]; valid=np.isfinite(x)&(x>=0)&(x<=1); good[var]=valid
                                        c=counts[var]; c['selected']+=len(x); c['valid']+=int(valid.sum())
                                        c['negative']+=int((np.isfinite(x)&(x<0)).sum()); c['nonfinite']+=int((~np.isfinite(x)).sum()); c['above_one']+=int((np.isfinite(x)&(x>1)).sum())
                                        c['sumw']+=float(w[valid].sum()); c['sumw2']+=float((w[valid]**2).sum())
                                        key=f'{cat}_{region}_{var}__{proc}'
                                        archive[key][0]+=np.histogram(x[valid],EDGES,weights=w[valid])[0]
                                        archive[key][1]+=np.histogram(x[valid],EDGES,weights=w[valid]**2)[0]
                                        if proc=='DY':
                                            off_counts[var]['sumw']+=float(wo[valid].sum()); off_counts[var]['sumw2']+=float((wo[valid]**2).sum())
                                            archive[f'{cat}_{region}_{var}__DY_no_pnet'][0]+=np.histogram(x[valid],EDGES,weights=wo[valid])[0]
                                            archive[f'{cat}_{region}_{var}__DY_no_pnet'][1]+=np.histogram(x[valid],EDGES,weights=wo[valid]**2)[0]
                            for var,c in counts.items():
                                assert c['selected']==c['valid']+c['negative']+c['nonfinite']+c['above_one']
                                rows.append(dict(era=era,category=cat,region=region,tag=tag,channel=channel,process=proc,variable=var,**c))
                            if proc=='DY':
                                for var,c in off_counts.items():
                                    row=counts[var].copy(); row.update(c)
                                    rows.append(dict(era=era,category=cat,region=region,tag=tag,channel=channel,process='DY_no_pnet',variable=var,**row))
            print(f'Loaded {cat} {region}',flush=True)
    assert len(sources)==168 and len(rows)==432
    for cat in ('merged',):
        for region in ('SR','CR'):
            for var in VARIABLES:
                for proc in (MC+('DY_no_pnet',) if region=='SR' else MC+('Data','DY_no_pnet')):
                    matching=[r for r in rows if (r['category'],r['region'],r['variable'],r['process'])==(cat,region,var,proc)]
                    np.testing.assert_allclose(archive[f'{cat}_{region}_{var}__{proc}'].sum(axis=1),[sum(r[k] for r in matching) for k in ('sumw','sumw2')],rtol=1e-10,atol=1e-6)
    assert not any('__Data' in k and '_SR_' in k for k in archive)
    np.savez_compressed(out/'histograms.npz',**archive)
    with (out/'validity_audit.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    (out/'metadata.json').write_text(json.dumps(dict(input=str(source),source_files=sources,validation_report_sha256=hashes,SR_Data_read=False,luminosity_display=138,
        validity='finite and 0<=score<=1; invalid scores excluded',SF_comparison='only DY toggled using weight and weight_no_pnet; other MC unchanged',
        selections='original era/category cuts, m2l2q>=300, 40<=mass2q<180; full CR includes bridge and VR',weighting='full nominal weight including DY ParticleNet SF; original era luminosities retained'),indent=2))
    print('PASS: all valid-event histogram yields and variances; SR Data absent',flush=True)
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
    dy=archive[f'{cat}_{region}_{var}__DY']; dyoff=archive[f'{cat}_{region}_{var}__DY_no_pnet']
    totaloff=total-dy[0]+dyoff[0]; erroff=np.sqrt(variance-dy[1]+dyoff[1])
    offlabel='Total, DY SF off'; onlabel='Total, DY SF on'
    ax.stairs(total,EDGES,color='#245b36',lw=2,label=onlabel)
    ax.stairs(totaloff,EDGES,color='#b32936',lw=2.3,ls='--',label=offlabel)
    maximum=float(max(np.max(total+err),np.max(totaloff+erroff))); centers=(EDGES[:-1]+EDGES[1:])/2
    if region=='CR':
        data=archive[f'{cat}_{region}_{var}__Data'][0]; de=poisson(data)
        ax.errorbar(centers,data,yerr=de,fmt='o',color='black',ms=5,lw=1.2,label='Data',zorder=5)
        maximum=max(maximum,float((data+de[1]).max())); good=total>0
        rel=np.divide(err,total,out=np.full_like(err,np.nan),where=good)
        ratio.fill_between(EDGES,step(1-rel),step(1+rel),step='post',color='#4c566a',alpha=.3,lw=0)
        ratio.errorbar(centers[good],data[good]/total[good],yerr=de[:,good]/total[good],fmt='o',color='black',ms=5,lw=1.2)
        go=totaloff>0
        ratio.errorbar(centers[go],data[go]/totaloff[go],yerr=de[:,go]/totaloff[go],fmt='s',mfc='none',color='#b32936',ms=5,lw=.9)
        relo=np.divide(erroff,totaloff,out=np.full_like(erroff,np.nan),where=go)
        ratio.stairs(1+relo,EDGES,color='#b32936',lw=.8,ls='--'); ratio.stairs(1-relo,EDGES,color='#b32936',lw=.8,ls='--')
        off=good&(data>2.5*total); ratio.plot(centers[off],np.full(int(off.sum()),2.43),'^',color='black',ms=7)
        ratio.axhline(1,color='.4',lw=1); ratio.set_ylim(0,2.5); ratio.set_ylabel('Data / MC',fontsize=27)
    positive=total[total>0]
    ax.set_yscale('log'); ax.set_ylim(.1,max(maximum,1.)*10000)
    hep.cms.label('Preliminary',data=region=='CR',lumi=138,com=13,ax=ax,fontsize=27)
    ax.text(.035,.95,f'Run 2, {cat.capitalize()}\nAll tags, '+r'$ee+\mu\mu$'+'\n'+('SR (blinded)' if region=='SR' else 'CR'),transform=ax.transAxes,va='top',fontsize=23)
    handles,labels=ax.get_legend_handles_labels(); lookup=dict(zip(labels,handles)); order=[LABELS[p] for p in MC]+[offlabel,band]+(['Data'] if region=='CR' else [])
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
    with PdfPages(out/'all_Run2_PNet_scores.pdf') as bundle:
        for region in ('SR','CR'):
            for cat in ('merged',):
                fig=plt.figure(figsize=(22,10)); outer=fig.add_gridspec(1,2,left=.06,right=.985,bottom=.16,top=.9,wspace=.32)
                for i,var in enumerate(VARIABLES):
                    one=draw(archive,cat,region,var)
                    for ext in ('png','pdf'): one.savefig(out/f'{cat}_{region}_{var}.{ext}',dpi=140)
                    bundle.savefig(one); plt.close(one)
                    if region=='CR':
                        sub=outer[i].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05); ax=fig.add_subplot(sub[0]); ratio=fig.add_subplot(sub[1],sharex=ax)
                    else: ax=fig.add_subplot(outer[i]); ratio=None
                    draw(archive,cat,region,var,(ax,ratio))
                for ext in ('png','pdf'): fig.savefig(out/f'{cat}_{region}_PNet_scores.{ext}',dpi=120)
                plt.close(fig)
    print('Completed 4 individual comparisons and 2 overview pages',flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--input',type=Path,required=True); p.add_argument('--output',type=Path,required=True); p.add_argument('--histograms',type=Path)
    a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
    if a.histograms:
        with np.load(a.histograms) as old: archive={k:old[k].copy() for k in old.files}
    else: archive=load(a.input,a.output)
    render(archive,a.output)
if __name__=='__main__': main()
