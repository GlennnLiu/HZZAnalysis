"""Data versus alpha*(Data-others)_CR+others_target, DY PNet off/on."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import uproot
from scipy.stats import chi2
from matplotlib.backends.backend_pdf import PdfPages
from plot_helpers import draw,plt,hep,TAGS


def primitive(source,process,category,tag,channel,region,target,edges):
    h=np.zeros((2,len(edges)-1))
    with uproot.open(source/f'{process}_{channel}_{category}_{tag}_2018_{region}.root') as f:
        for a in f['tree'].iterate(['m2l2q','mass2q','weight'],step_size=100000,library='np'):
            q=a['mass2q']; mask=(a['m2l2q']>=300)&(q>=40)&(q<180)
            if region=='CR':
                assert ((q<70)|(q>105)).all()
                if target=='VR': mask&=~((q>135)&(q<165))
            else: assert (((q>70)&(q<105)) if target=='SR' else ((q>135)&(q<165))).all()
            w=a['weight'][mask]
            if process=='Data': np.testing.assert_array_equal(w,1)
            for i,v in enumerate((w,w*w)): h[i]+=np.histogram(a['m2l2q'][mask],bins=edges,weights=v)[0]
    return h


def prediction(alpha,z,other):
    p=alpha[:2]*z[0]+other[0]
    var=alpha[:2]**2*z[1]+z[0]**2*alpha[2:4]+other[1]
    cov=alpha[0]*alpha[1]*z[1]+z[0]**2*alpha[4]+other[1]
    return np.vstack((p,var,cov))


def statistics(pred,data,mask):
    rows=[]; valid=np.isfinite(pred).all(axis=0)
    for i,method in enumerate(('no_pnet','with_pnet')):
        v=data[1]+pred[i+2]
        value=float(np.sum((data[0,mask]-pred[i,mask])**2/v[mask])); n=int(mask.sum())
        d=float(data[0,valid].sum()); p=float(pred[i,valid].sum()); pv=float(pred[i+2,valid].sum())
        rows.append(dict(method=method,chi2=value,ndof=n,pvalue=float(chi2.sf(value,n)) if n else None,
            data=d,prediction=p,prediction_stat_error=np.sqrt(pv),data_over_bkg=d/p if p>0 else None,
            yield_pull=(d-p)/np.sqrt(d+pv),valid_bins=int(valid.sum()),unavailable_data_yield=float(data[0,~valid].sum())))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--input',type=Path,required=True); p.add_argument('--alpha-archive',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    args=p.parse_args(); args.output.mkdir(parents=True,exist_ok=True); plt.style.use(hep.style.CMS)
    cache={}; archive={}; summary=[]
    with np.load(args.alpha_archive) as old:
        for region in ('SR','VR'):
            for category in ('resolved','merged'):
                for tag in TAGS:
                    for channel in ('ee','mumu'):
                        key=f'{category}_{tag}_{channel}_{region}'; edges=old[key+'__edges']; alpha=old[key+'__alpha']
                        source={process:primitive(args.input,process,category,tag,channel,'CR',region,edges) for process in ('Data','TTplusWW','VZ')}
                        target={process:primitive(args.input,process,category,tag,channel,region,region,edges) for process in ('Data','TTplusWW','VZ')}
                        z=np.stack((source['Data'][0]-source['TTplusWW'][0]-source['VZ'][0],sum(h[1] for h in source.values())))
                        other=target['TTplusWW']+target['VZ']; data=target['Data']
                        pred=prediction(alpha,z,other)
                        valid=np.isfinite(pred).all(axis=0)
                        mask=old[key+'__stable']&(data[0]>=5)&valid&(pred[0]>0)&(pred[1]>0)
                        assert np.all(pred[2:4,valid]>=0)
                        assert np.all(pred[4,valid]**2<=pred[2,valid]*pred[3,valid]+1e-6)
                        if category=='resolved': np.testing.assert_array_equal(pred[0,valid],pred[1,valid])
                        rows=statistics(pred,data,mask)
                        for row in rows: summary.append(dict(category=category,tag=tag,channel=channel,region=region,**row))
                        zeros=np.zeros_like(pred)
                        parts=dict(DY=pred,TTplusWW=zeros,VZ=zeros,Data=np.tile(data[0],(5,1)))
                        cache[key]=edges,parts
                        for field,values in (('edges',edges),('alpha',alpha),('z_source',z),('other_target',other),('data',data),('prediction',pred),('chi2_mask',mask)):
                            archive[key+'__'+field]=values
                        for process in source: archive[key+'__source_'+process]=source[process]; archive[key+'__target_'+process]=target[process]
            print(f'Computed {region}',flush=True)
    # Combined diagnostic chi2 is the sum of separate-channel chi2s, not a pooled alpha.
    for region in ('SR','VR'):
        for category in ('resolved','merged'):
            for tag in TAGS:
                for method in ('no_pnet','with_pnet'):
                    rows=[r for r in summary if r['region']==region and r['category']==category and r['tag']==tag and r['method']==method and r['channel'] in ('ee','mumu')]
                    d=sum(r['data'] for r in rows); pred=sum(r['prediction'] for r in rows)
                    value=sum(r['chi2'] for r in rows); n=sum(r['ndof'] for r in rows)
                    summary.append(dict(category=category,tag=tag,channel='combined_channel_tests',region=region,method=method,chi2=value,ndof=n,pvalue=float(chi2.sf(value,n)) if n else None,
                        data=d,prediction=pred,prediction_stat_error=np.sqrt(sum(r['prediction_stat_error']**2 for r in rows)),data_over_bkg=d/pred if pred>0 else None,
                        yield_pull=(d-pred)/np.sqrt(d+sum(r['prediction_stat_error']**2 for r in rows)),valid_bins=sum(r['valid_bins'] for r in rows),unavailable_data_yield=sum(r['unavailable_data_yield'] for r in rows)))
    with (args.output/'compatibility.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0])); w.writeheader(); w.writerows(summary)
    np.savez_compressed(args.output/'histograms.npz',**archive)
    (args.output/'metadata.json').write_text(json.dumps(dict(input=str(args.input),alpha_archive=str(args.alpha_archive),alpha_archive_sha256=hashlib.sha256(args.alpha_archive.read_bytes()).hexdigest(),statistics='First-order stat only; Data>=5 and DY Neff>=5 in source/target for both choices; common masks; no fitted normalization'),indent=2))
    for r in summary:
        if r['channel']=='combined_channel_tests': print(json.dumps(r),flush=True)
    for region in ('SR','VR'):
        with PdfPages(args.output/f'all_{region}_alpha_DataBkg_PNet.pdf') as bundle:
            for channel in ('ee','mumu'):
                fig=plt.figure(figsize=(33,20)); outer=fig.add_gridspec(2,3,left=.065,right=.99,bottom=.09,top=.95,hspace=.40,wspace=.35)
                for row,category in enumerate(('resolved','merged')):
                    for col,tag in enumerate(TAGS):
                        key=f'{category}_{tag}_{channel}_{region}'; edges,parts=cache[key]
                        kwargs=dict(edges_override=edges,channel=channel,region=region)
                        one=draw(parts,category,tag,'m2l2q',**kwargs)
                        for ext in ('png','pdf'): one.savefig(args.output/f'{key}_alpha_DataBkg_PNet.{ext}',dpi=140)
                        bundle.savefig(one); plt.close(one)
                        sub=outer[row,col].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05)
                        ax=fig.add_subplot(sub[0]); lower=fig.add_subplot(sub[1],sharex=ax)
                        draw(parts,category,tag,'m2l2q',(ax,lower),**kwargs)
                for ext in ('png','pdf'): fig.savefig(args.output/f'six_categories_{channel}_{region}_alpha_DataBkg_PNet.{ext}',dpi=120)
                plt.close(fig)


if __name__=='__main__': main()
