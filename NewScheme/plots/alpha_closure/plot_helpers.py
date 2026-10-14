"""2018 SR observed Data versus full MC, varying DY ParticleNet only."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
import argparse
import csv
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

MASS=np.array([300,400,500,600,700,800,900,1000,1200,1400,1600,2000,2500,3000,np.inf])
TAGS=('untag','btag','vbftag')
LABELS={'untag':'untagged','btag':'b-tag','vbftag':'VBF-tag'}


def poisson_errors(n):
    low=np.zeros_like(n); positive=n>0
    low[positive]=.5*chi2.ppf(.15865525393145707,2*n[positive])
    high=.5*chi2.ppf(.8413447460685429,2*(n+1))
    return np.stack((n-low,high-n))


def load(source,category,tag):
    hist={kind:{p:np.zeros((5,len(edges)-1)) for p in ('DY','TTplusWW','VZ','Data')}
          for kind,edges in (('mqq',np.arange(70.,110.,5.)),('m2l2q',MASS))}
    for process in ('DY','TTplusWW','VZ','Data'):
        for channel in ('ee','mumu'):
            path=source/f'{process}_{channel}_{category}_{tag}_2018_SR.root'
            with uproot.open(path) as handle:
                for a in handle['tree'].iterate(['mass2q','m2l2q','weight','weight_no_pnet','pnet_sf'],step_size=100000,library='np'):
                    mask=a['m2l2q']>=300
                    assert ((a['mass2q']>70)&(a['mass2q']<105)).all()
                    w=a['weight'][mask]
                    u=a['weight_no_pnet'][mask] if process=='DY' else w
                    if process=='Data': np.testing.assert_array_equal(w,np.ones(len(w)))
                    if category=='resolved' and process=='DY': np.testing.assert_array_equal(w,u)
                    for kind,edges in (('mqq',np.arange(70.,110.,5.)),('m2l2q',MASS)):
                        x=a['mass2q' if kind=='mqq' else 'm2l2q'][mask]
                        for i,weights in enumerate((u,w,u*u,w*w,u*w)):
                            hist[kind][process][i]+=np.histogram(x,bins=edges,weights=weights)[0]
    return hist


def calculate(parts):
    mc=sum(parts[p] for p in ('DY','TTplusWW','VZ'))
    d=parts['Data'][0]
    np.testing.assert_array_equal(d,parts['Data'][2])
    neff0=np.divide(mc[0]**2,mc[2],out=np.zeros_like(d),where=mc[2]>0)
    neff1=np.divide(mc[1]**2,mc[3],out=np.zeros_like(d),where=mc[3]>0)
    mask=(d>=5)&(mc[0]>0)&(mc[1]>0)&(neff0>=5)&(neff1>=5)
    rows=[]
    for i,label in enumerate(('DY_no_pnet','DY_with_pnet')):
        variance=d+mc[i+2]
        value=float(np.sum((d[mask]-mc[i,mask])**2/variance[mask]))
        n=int(mask.sum()); prediction=float(mc[i].sum())
        rows.append(dict(method=label,data=float(d.sum()),mc=prediction,dy=float(parts['DY'][i].sum()),
            other_mc=float((parts['TTplusWW'][i]+parts['VZ'][i]).sum()),data_over_mc=float(d.sum()/prediction),
            yield_pull=float((d.sum()-prediction)/np.sqrt(d.sum()+mc[i+2].sum())),
            chi2=value,ndof=n,pvalue=float(chi2.sf(value,n)) if n else None))
    return mc,d,mask,rows


def draw(parts,category,tag,kind,axes=None,edges_override=None,channel='combined',region='SR'):
    mc,d,mask,rows=calculate(parts)
    if axes is None:
        fig,(ax,lower)=plt.subplots(2,1,figsize=(11,10),sharex=True,gridspec_kw={'height_ratios':[3.3,1],'hspace':.05})
        fig.subplots_adjust(left=.18,right=.97,top=.90,bottom=.16)
    else: ax,lower=axes; fig=ax.figure
    edges=edges_override.copy() if edges_override is not None else (np.arange(70.,110.,5.) if kind=='mqq' else MASS.copy())
    if np.isinf(edges[-1]): edges[-1]=4000
    scale=np.ones(len(d)) if kind=='mqq' else 100/np.diff(edges)
    centers=(edges[:-1]+edges[1:])/2
    data_err=poisson_errors(d)
    ratio_upper=[2.5]; ratio_lower=[0.]
    for i,(color,label) in enumerate((('#2765a3','Bkg.: alpha without PNet SF'),('#d1493f','Bkg.: alpha with PNet SF'))):
        y=mc[i]*scale; err=np.sqrt(mc[i+2])*scale
        ax.stairs(y,edges,color=color,lw=2.5,label=label)
        ax.fill_between(edges,np.r_[y-err,(y-err)[-1]],np.r_[y+err,(y+err)[-1]],step='post',alpha=.14,color=color,lw=0)
        good=mc[i]>0
        lower.errorbar(centers[good],d[good]/mc[i,good],yerr=data_err[:,good]/mc[i,good],fmt='o',color=color,ms=5)
        off=good&(d>2.5*mc[i])
        lower.plot(centers[off],np.full(int(off.sum()),2.43),'^',color=color,ms=7)
        ratio_upper.extend(((d[good]+data_err[1,good])/mc[i,good]).tolist())
        relative=np.divide(np.sqrt(mc[i+2]),mc[i],out=np.full_like(d,np.nan),where=good)
        ratio_upper.extend((1+relative[good]).tolist()); ratio_lower.extend((1-relative[good]).tolist())
        lower.fill_between(edges,np.r_[1-relative,(1-relative)[-1]],np.r_[1+relative,(1+relative)[-1]],step='post',color=color,alpha=.12,lw=0)
    ax.errorbar(centers,d*scale,yerr=data_err*scale,fmt='o',color='black',ms=5,label='Observed data')
    plotted=np.r_[d*scale,mc[:2].ravel()*np.tile(scale,2)]; positive=plotted[plotted>0]
    ax.set_yscale('log'); ax.set_ylim(max(.001,float(positive.min())*.2),max(1.,float(np.nanmax(plotted)))*80)
    hep.cms.label('Preliminary',data=True,lumi=59.83,com=13,ax=ax,fontsize=27)
    ax.text(.035,.95,f'2018, {category.capitalize()}\n{LABELS[tag]}, {region}\n'+(r'$ee$' if channel=='ee' else r'$\mu\mu$'),transform=ax.transAxes,va='top',fontsize=23)
    ax.legend(loc='upper right',fontsize=20,frameon=False)
    ax.set_ylabel('Events / 5 GeV' if kind=='mqq' else 'Events / 100 GeV',fontsize=32)
    lower.axhline(1,color='.4',lw=1); lower.set_ylim(0,2.5)
    lower.set_ylabel('Data / Bkg.',fontsize=27)
    lower.set_xlabel((r'$m_{jj}$ [GeV]' if category=='resolved' else r'$m_J$ [GeV]') if kind=='mqq' else r'$m_{2\ell2q}$ [GeV]',fontsize=34)
    lower.set_xlim(edges[0],edges[-1]); lower.set_xticks([x for x in (300,1000,2000,3000,3500,4000) if edges[0]<=x<=edges[-1]])
    for axis in (ax,lower): axis.tick_params(direction='in',top=True,right=True,labelsize=23); axis.grid(False)
    ax.tick_params(labelbottom=False)
    return fig


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--input',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    p.add_argument('--histograms',type=Path,help='Reuse saved histograms without rereading trees')
    args=p.parse_args(); args.output.mkdir(parents=True,exist_ok=True); plt.style.use(hep.style.CMS)
    cache={}; archive={}; summary=[]
    for category in ('resolved','merged'):
        report=json.loads((args.input/f'validation_2018_{category}.json').read_text()); assert report['status']=='passed' and not report['pilot']
        for tag in TAGS:
            if args.histograms:
                with np.load(args.histograms) as old:
                    h={kind:{process:old[f'{category}_{tag}_{kind}__{process}'].copy()
                             for process in ('DY','TTplusWW','VZ','Data')} for kind in ('mqq','m2l2q')}
            else: h=load(args.input,category,tag)
            np.testing.assert_allclose(h['mqq']['Data'].sum(axis=1),h['m2l2q']['Data'].sum(axis=1))
            for kind,parts in h.items():
                key=f'{category}_{tag}_{kind}'; cache[key]=parts
                for process,values in parts.items(): archive[key+'__'+process]=values
                mc,d,mask,rows=calculate(parts); archive[key+'__chi2_mask']=mask
                for row in rows: summary.append(dict(category=category,tag=tag,variable=kind,**row))
    with PdfPages(args.output/'all_SR_DataMC_PNet.pdf') as bundle:
        for kind in ('mqq','m2l2q'):
            fig=plt.figure(figsize=(33,20)); outer=fig.add_gridspec(2,3,left=.065,right=.99,bottom=.09,top=.95,hspace=.40,wspace=.35)
            for row,category in enumerate(('resolved','merged')):
                for col,tag in enumerate(TAGS):
                    key=f'{category}_{tag}_{kind}'; parts=cache[key]
                    one=draw(parts,category,tag,kind)
                    for ext in ('png','pdf'): one.savefig(args.output/f'{key}_SR_DataMC_PNet.{ext}',dpi=140)
                    bundle.savefig(one); plt.close(one)
                    sub=outer[row,col].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05)
                    ax=fig.add_subplot(sub[0]); lower=fig.add_subplot(sub[1],sharex=ax)
                    draw(parts,category,tag,kind,(ax,lower))
            for ext in ('png','pdf'): fig.savefig(args.output/f'six_categories_SR_{kind}_DataMC_PNet.{ext}',dpi=120)
            plt.close(fig)
    with (args.output/'compatibility.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(summary[0])); writer.writeheader(); writer.writerows(summary)
    np.savez_compressed(args.output/'histograms.npz',**archive,m2l2q_edges=MASS,mqq_edges=np.arange(70.,110.,5.))
    (args.output/'results.json').write_text(json.dumps(summary,indent=2))
    for r in summary:
        if r['variable']=='m2l2q': print(json.dumps(r),flush=True)


if __name__=='__main__': main()
