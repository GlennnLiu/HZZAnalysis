"""Paired DY ParticleNet effect on SR/full-CR and VR/held-out-CR alpha."""
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
import uproot

TAGS=('untag','btag','vbftag')
LABELS={'untag':'untagged','btag':'b-tag','vbftag':'VBF-tag'}


def paired_alpha(t,c):
    """Rows are noSF, SF, var(noSF), var(SF), paired covariance."""
    valid=(t[0]>0)&(t[1]>0)&(c[0]>0)&(c[1]>0)
    out=np.full((8,t.shape[1]),np.nan)
    u,w,vu,vw,ct=t[:,valid]; b,d,vb,vd,cc=c[:,valid]
    a0=u/b; a1=w/d
    v0=vu/b**2+u*u*vb/b**4; v1=vw/d**2+w*w*vd/d**4
    cov=ct/(b*d)+u*w*cc/(b*b*d*d)
    r=a1/a0
    vr=r*r*(vw/w**2+vu/u**2-2*ct/(u*w)+vd/d**2+vb/b**2-2*cc/(b*d))
    delta_var=v0+v1-2*cov
    assert np.all(vr>=-1e-8*np.maximum(1.,r*r))
    assert np.all(delta_var>=-1e-8*np.maximum(1.,v0+v1))
    out[:,valid]=np.stack((a0,a1,v0,v1,cov,r,np.maximum(vr,0),np.maximum(delta_var,0)))
    return out,valid


def load(source,category,tag,channel,target,edges):
    result=[]; totals=[]
    for region in (target,'CR'):
        h=np.zeros((5,len(edges)-1)); total=np.zeros((5,1))
        with uproot.open(source/f'DY_{channel}_{category}_{tag}_2018_{region}.root') as handle:
            for a in handle['tree'].iterate(['m2l2q','mass2q','weight','weight_no_pnet','pnet_sf'],step_size=100000,library='np'):
                q=a['mass2q']; mask=(a['m2l2q']>=300)&(q>=40)&(q<180)
                if region=='CR':
                    assert ((q<70)|(q>105)).all()
                    if target=='VR': mask&=~((q>135)&(q<165))
                else: assert (((q>70)&(q<105)) if target=='SR' else ((q>135)&(q<165))).all()
                w=a['weight'][mask]; u=a['weight_no_pnet'][mask]
                np.testing.assert_allclose(w,u*a['pnet_sf'][mask],rtol=1e-14)
                if category=='resolved': np.testing.assert_array_equal(w,u)
                for i,values in enumerate((u,w,u*u,w*w,u*w)):
                    h[i]+=np.histogram(a['m2l2q'][mask],bins=edges,weights=values)[0]
                    total[i,0]+=values.sum()
        result.append(h); totals.append(total)
    return result,totals


def stable_mask(t,c,valid):
    mask=valid.copy()
    for h in (t,c):
        for i in (0,1):
            mask&=np.divide(h[i]**2,h[i+2],out=np.zeros_like(h[i]),where=h[i+2]>0)>=5
    return mask


def draw(edges,values,category,tag,channel,region,axes=None):
    if axes is None:
        fig,(ax,lower)=plt.subplots(2,1,figsize=(11,10),sharex=True,gridspec_kw={'height_ratios':[3.3,1],'hspace':.05})
        fig.subplots_adjust(left=.18,right=.97,top=.90,bottom=.16)
    else: ax,lower=axes; fig=ax.figure
    centers=(edges[:-1]+edges[1:])/2; valid=np.isfinite(values[0])
    for i,(color,label) in enumerate((('#2765a3','Without ParticleNet SF'),('#d1493f','With ParticleNet SF'))):
        ax.errorbar(centers[valid],values[i,valid],xerr=np.diff(edges)[valid]/2,yerr=np.sqrt(values[i+2,valid]),fmt='o',color=color,ms=5,lw=1.2,label=label)
    err=np.sqrt(values[6]); lower.axhline(1,color='.4',lw=1)
    lower.errorbar(centers[valid],values[5,valid],yerr=err[valid],fmt='o',color='#d1493f',ms=5)
    lower.set_ylim(min(.9,float(np.nanmin(values[5]-err))-.05),max(1.1,float(np.nanmax(values[5]+err))+.05))
    lo=min(0.,float(np.nanmin(values[:2]-np.sqrt(values[2:4])))); hi=float(np.nanmax(values[:2]+np.sqrt(values[2:4])))
    ax.set_ylim(lo-.05*max(hi-lo,.1),hi+.65*max(hi-lo,.1))
    hep.cms.label('Preliminary',data=False,lumi=59.83,com=13,ax=ax,fontsize=27)
    ax.text(.035,.95,f'2018, {category.capitalize()}\n{LABELS[tag]}, '+(r'$ee$' if channel=='ee' else r'$\mu\mu$')+f'\n{region} / '+('full CR' if region=='SR' else 'CR excluding VR'),transform=ax.transAxes,va='top',fontsize=23)
    ax.legend(loc='upper right',fontsize=20,frameon=False)
    ax.set_ylabel(r'$\alpha$',fontsize=34); lower.set_ylabel('With / without',fontsize=25)
    lower.set_xlabel(r'$m_{2\ell2q}$ [GeV]',fontsize=34); lower.set_xlim(edges[0],edges[-1])
    lower.set_xticks([x for x in (300,1000,2000,3000,3500,4000) if edges[0]<=x<=edges[-1]])
    for axis in (ax,lower): axis.tick_params(direction='in',top=True,right=True,labelsize=23); axis.grid(False)
    ax.tick_params(labelbottom=False)
    return fig


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--input',type=Path,required=True); p.add_argument('--alpha-dir',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    args=p.parse_args(); args.output.mkdir(parents=True,exist_ok=True); plt.style.use(hep.style.CMS)
    archive={}; cache={}; summary=[]; bins=[]; audits=[]; total_cache={}
    for region in ('SR','VR'):
        path=args.alpha_dir/f'AlphaHisto_{region}_2018.root'
        audits.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        with uproot.open(path) as old:
            for category in ('resolved','merged'):
                report=json.loads((args.input/f'validation_2018_{category}.json').read_text()); assert report['status']=='passed' and not report['pilot']
                for tag in TAGS:
                    for channel in ('ee','mumu'):
                        key=f'{category}_{tag}_{channel}_{region}'
                        rootkey=f'{category if category=="resolved" else "merged_tag"}_{"isEE" if channel=="ee" else "isMuMu"}_{tag}'
                        raw=old[rootkey].axis().edges()
                        edges=np.r_[300.,raw[raw>300]]
                        (t,c),totals=load(args.input,category,tag,channel,region,edges)
                        values,valid=paired_alpha(t,c); stable=stable_mask(t,c,valid)
                        if category=='resolved':
                            np.testing.assert_array_equal(values[0,valid],values[1,valid]); np.testing.assert_allclose(values[5,valid],1)
                        cache[key]=edges,values
                        total_cache[category,tag,channel,region]=totals
                        for field,array in (('edges',edges),('target',t),('source',c),('alpha',values),('valid',valid),('stable',stable)): archive[key+'__'+field]=array
                        archive[key+'__target_total']=totals[0]; archive[key+'__source_total']=totals[1]
                        integrated,_=paired_alpha(*totals)
                        r=integrated[5,0]
                        summary.append(dict(category=category,tag=tag,channel=channel,region=region,
                            alpha_noSF=integrated[0,0],alpha_SF=integrated[1,0],relative_change_percent=100*(r-1),ratio_stat_error=np.sqrt(integrated[6,0]),
                            stable_bins=int(stable.sum()),max_abs_bin_change_percent=float(np.max(np.abs(values[5,stable]-1))*100) if stable.any() else None))
                        for i in range(len(edges)-1):
                            bins.append(dict(key=key,low=edges[i],high=edges[i+1],valid=bool(valid[i]),stable=bool(stable[i]),alpha_noSF=values[0,i],alpha_SF=values[1,i],ratio=values[5,i],ratio_stat_error=np.sqrt(values[6,i])))
    for region in ('SR','VR'):
        for category in ('resolved','merged'):
            for tag in TAGS:
                totals=[sum(total_cache[category,tag,ch,region][i] for ch in ('ee','mumu')) for i in (0,1)]
                values,_=paired_alpha(*totals)
                row=dict(category=category,tag=tag,channel='combined_integrated_DY_ratio',region=region,alpha_noSF=values[0,0],alpha_SF=values[1,0],relative_change_percent=100*(values[5,0]-1),ratio_stat_error=np.sqrt(values[6,0]),stable_bins=None,max_abs_bin_change_percent=None)
                summary.append(row)
                if category=='merged': print(json.dumps(row),flush=True)
        with PdfPages(args.output/f'all_{region}_alpha_PNet.pdf') as bundle:
            for channel in ('ee','mumu'):
                fig=plt.figure(figsize=(33,20)); outer=fig.add_gridspec(2,3,left=.065,right=.99,bottom=.09,top=.95,hspace=.40,wspace=.35)
                for row,category in enumerate(('resolved','merged')):
                    for col,tag in enumerate(TAGS):
                        key=f'{category}_{tag}_{channel}_{region}'; edges,values=cache[key]
                        one=draw(edges,values,category,tag,channel,region)
                        for ext in ('png','pdf'): one.savefig(args.output/f'{key}_alpha_PNet.{ext}',dpi=140)
                        bundle.savefig(one); plt.close(one)
                        sub=outer[row,col].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05)
                        ax=fig.add_subplot(sub[0]); lower=fig.add_subplot(sub[1],sharex=ax)
                        draw(edges,values,category,tag,channel,region,(ax,lower))
                for ext in ('png','pdf'): fig.savefig(args.output/f'six_categories_{channel}_{region}_alpha_PNet.{ext}',dpi=120)
                plt.close(fig)
    for filename,rows in (('summary.csv',summary),('bins.csv',bins)):
        with (args.output/filename).open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    np.savez_compressed(args.output/'histograms.npz',**archive)
    (args.output/'input_audit.json').write_text(json.dumps(audits,indent=2))


if __name__=='__main__': main()
