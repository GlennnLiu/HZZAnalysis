"""All-era/category/channel direct MC Z kinematics; SR Data never read."""
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
COLORS={'VZ':'#f4b183','TTplusWW':'#b1a2c5','DY':'#80c17d'}
LABELS={'VZ':'VZ','TTplusWW':r'$t\bar{t}$ + WW + tZq','DY':'DY + jets'}
PT=np.r_[np.arange(0.,1050.,50.),1200.,1500.,2000.,np.inf]


def edges_for(variable,region):
    if variable.startswith('pt'): return PT.copy()
    if variable=='mass2l': return np.arange(40.,142.,2.)
    return np.arange(70.,110.,5.) if region=='SR' else np.arange(40.,185.,5.)


def visible_bins(variable,region,edges):
    if variable=='mass2q' and region=='CR': return (edges[1:]<=70)|(edges[:-1]>=105)
    return np.ones(len(edges)-1,dtype=bool)


def poisson_errors(n):
    low=np.zeros_like(n); good=n>0
    low[good]=.5*chi2.ppf(.15865525393145707,2*n[good])
    high=.5*chi2.ppf(.8413447460685429,2*(n+1))
    return np.stack((n-low,high-n))


def load(source):
    archive={}; audits=[]; totals=[]; report_hashes={}
    for category in ('resolved','merged'):
        for era in ERAS:
            path=source/f'validation_{era}_{category}.json'; r=json.loads(path.read_text())
            assert r['status']=='passed' and not r['pilot'] and (source/f'{era}_{category}.done').is_file()
            report_hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        for region in ('SR','CR'):
            processes=MC if region=='SR' else MC+('Data',)
            variables=('pt2l','mass2l','pt2q','mass2q')
            for variable in variables:
                edges=edges_for(variable,region); prefix=f'{category}_{region}_{variable}'
                archive[prefix+'__edges']=edges
                for process in processes: archive[prefix+'__'+process]=np.zeros((2,len(edges)-1))
            for era in ERAS:
                for tag in TAGS:
                    for channel in ('ee','mumu'):
                        for process in processes:
                            path=source/f'{process}_{channel}_{category}_{tag}_{era}_{region}.root'
                            audits.append(dict(file=str(path),size=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
                            n=0; sumw=0.; sumw2=0.; bridge=0; vr=0
                            with uproot.open(path) as handle:
                                for a in handle['tree'].iterate(list(variables)+['m2l2q','weight'],step_size=100000,library='np'):
                                    q=a['mass2q']; mask=(a['m2l2q']>=300)&(q>=40)&(q<180)
                                    assert (((q>70)&(q<105)) if region=='SR' else ((q<70)|(q>105))).all()
                                    w=a['weight'][mask]
                                    if process=='Data': np.testing.assert_array_equal(w,np.ones(len(w)))
                                    n+=int(mask.sum()); sumw+=float(w.sum()); sumw2+=float((w*w).sum())
                                    bridge+=int((mask&(q>105)&(q<135)).sum()); vr+=int((mask&(q>135)&(q<165)).sum())
                                    for variable in variables:
                                        prefix=f'{category}_{region}_{variable}'; x=a[variable][mask]; edges=archive[prefix+'__edges']
                                        assert np.isfinite(x).all() and np.isfinite(w).all()
                                        for i,values in enumerate((w,w*w)):
                                            archive[prefix+'__'+process][i]+=np.histogram(x,bins=edges,weights=values)[0]
                            totals.append(dict(era=era,category=category,region=region,tag=tag,channel=channel,process=process,entries=n,sumw=sumw,sumw2=sumw2,bridge_entries=bridge,VR_entries=vr))
            print(f'Loaded {category} {region}',flush=True)
    assert len(audits)==336 and not any('_SR.root' in r['file'] and '/Data_' in r['file'] for r in audits)
    return archive,audits,totals,report_hashes


def draw(archive,category,region,variable,axes=None):
    key=f'{category}_{region}_{variable}'; edges=archive[key+'__edges']; display=edges.copy()
    if np.isinf(display[-1]): display[-1]=2500.
    norm=50. if variable.startswith('pt') else (2. if variable=='mass2l' else 5.)
    scale=norm/np.diff(display); mask=visible_bins(variable,region,edges)
    parts={p:archive[key+'__'+p] for p in MC}
    total=sum(h[0] for h in parts.values()); variance=sum(h[1] for h in parts.values()); error=np.sqrt(variance)
    if axes is None:
        if region=='CR':
            fig,(ax,lower)=plt.subplots(2,1,figsize=(11,10),sharex=True,gridspec_kw={'height_ratios':[3.3,1],'hspace':.05})
            fig.subplots_adjust(left=.18,right=.97,top=.90,bottom=.16)
        else:
            fig,ax=plt.subplots(figsize=(11,9)); lower=None
            fig.subplots_adjust(left=.18,right=.97,top=.90,bottom=.16)
    else: ax,lower=axes; fig=ax.figure
    bottom=np.zeros_like(total)
    for process in MC:
        values=parts[process][0]*scale
        ax.bar(display[:-1][mask],values[mask],width=np.diff(display)[mask],bottom=bottom[mask],align='edge',color=COLORS[process],edgecolor='none',label=LABELS[process])
        bottom+=values
    step=lambda x:np.r_[x,x[-1]]
    low=np.where(mask,(total-error)*scale,np.nan); high=np.where(mask,(total+error)*scale,np.nan)
    ax.fill_between(display,step(low),step(high),step='post',color='#4c566a',alpha=.3,lw=0,label=r'Bkg. MC stat. ($1\sigma$)')
    maximum=float(np.max(high[mask])); positive=total[mask]*scale[mask]; positive=positive[positive>0]
    if region=='CR':
        data=archive[key+'__Data'][0]; err=poisson_errors(data); good=mask&(total>0)
        ax.errorbar((display[:-1]+display[1:])[mask]/2,data[mask]*scale[mask],yerr=err[:,mask]*scale[mask],fmt='o',color='black',ms=5,lw=1.2,label='Data',zorder=5)
        maximum=max(maximum,float(np.max((data[mask]+err[1,mask])*scale[mask])))
        relative=np.divide(error,total,out=np.full_like(total,np.nan),where=good)
        lower.fill_between(display,step(1-relative),step(1+relative),step='post',color='#4c566a',alpha=.3,lw=0)
        centers=(display[:-1]+display[1:])/2
        lower.errorbar(centers[good],data[good]/total[good],yerr=err[:,good]/total[good],fmt='o',color='black',ms=5,lw=1.2)
        off=good&(data>2.5*total); lower.plot(centers[off],np.full(int(off.sum()),2.43),'^',color='black',ms=7)
        lower.axhline(1,color='.4',lw=1); lower.set_ylim(0,2.5); lower.set_ylabel('Data / MC',fontsize=27)
    headroom=1000 if variable.startswith('pt') else 100
    ax.set_yscale('log'); ax.set_ylim(max(.001,float(positive.min())*.2),max(maximum,1.)*headroom)
    hep.cms.label('Preliminary',data=region=='CR',lumi=138,com=13,ax=ax,fontsize=27)
    ax.text(.035,.95,f'Run 2, {category.capitalize()}\nAll tags, '+r'$ee+\mu\mu$'+f'\n'+('SR (blinded)' if region=='SR' else 'CR'),transform=ax.transAxes,va='top',fontsize=23)
    handles,labels=ax.get_legend_handles_labels()
    order=[LABELS[p] for p in MC]+[r'Bkg. MC stat. ($1\sigma$)']+(['Data'] if region=='CR' else [])
    lookup=dict(zip(labels,handles))
    ax.legend([lookup[label] for label in order],order,loc='upper right',fontsize=22,frameon=False,handlelength=1.2,labelspacing=.3)
    ax.set_ylabel(f'Events / {norm:g} GeV',fontsize=32)
    q='jj' if category=='resolved' else 'J'
    label={'pt2l':r'$p_{T}^{\ell\ell}$ [GeV]','mass2l':r'$m_{\ell\ell}$ [GeV]',
           'pt2q':rf'$p_{{T}}^{{{q}}}$ [GeV]','mass2q':rf'$m_{{{q}}}$ [GeV]'}[variable]
    bottom_axis=lower if lower is not None else ax
    bottom_axis.set_xlabel(label,fontsize=34); bottom_axis.set_xlim(display[0],display[-1])
    if variable.startswith('pt'):
        bottom_axis.set_xticks([0,500,1000,1500,2000,2500])
        ax.text(.98,.47,r'Last bin: $\geq2000$ GeV',transform=ax.transAxes,ha='right',fontsize=17)
    elif variable=='mass2q': bottom_axis.set_xticks([70,80,90,100,105] if region=='SR' else [40,70,105,140,180])
    for axis in (ax,lower):
        if axis is not None: axis.tick_params(direction='in',top=True,right=True,labelsize=23); axis.grid(False)
    if lower is not None: ax.tick_params(labelbottom=False)
    return fig


def render(archive,out):
    plt.style.use(hep.style.CMS)
    with PdfPages(out/'all_Run2_Z_kinematics.pdf') as bundle:
        for region in ('SR','CR'):
            for category in ('resolved','merged'):
                fig=plt.figure(figsize=(22,20)); outer=fig.add_gridspec(2,2,left=.095,right=.985,bottom=.08,top=.95,hspace=.38,wspace=.35)
                for slot,variable in enumerate(('pt2l','mass2l','pt2q','mass2q')):
                    key=f'{category}_{region}_{variable}'
                    one=draw(archive,category,region,variable)
                    for ext in ('png','pdf'): one.savefig(out/f'{key}.{ext}',dpi=140)
                    bundle.savefig(one); plt.close(one)
                    if region=='CR':
                        sub=outer[slot//2,slot%2].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05)
                        ax=fig.add_subplot(sub[0]); lower=fig.add_subplot(sub[1],sharex=ax)
                    else: ax=fig.add_subplot(outer[slot//2,slot%2]); lower=None
                    draw(archive,category,region,variable,(ax,lower))
                for ext in ('png','pdf'): fig.savefig(out/f'{category}_{region}_four_variables.{ext}',dpi=120)
                plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--input',type=Path,required=True); p.add_argument('--output',type=Path,required=True); p.add_argument('--histograms',type=Path)
    args=p.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    if args.histograms:
        with np.load(args.histograms) as old: archive={k:old[k].copy() for k in old.files}
    else:
        archive,audit,totals,hashes=load(args.input)
        np.savez_compressed(args.output/'histograms.npz',**archive)
        with (args.output/'input_yields.csv').open('w',newline='') as handle:
            w=csv.DictWriter(handle,fieldnames=list(totals[0])); w.writeheader(); w.writerows(totals)
        (args.output/'metadata.json').write_text(json.dumps(dict(input=str(args.input),eras=list(ERAS),channels=['ee','mumu'],tags=list(TAGS),
            selection='m2l2q>=300, 40<=mass2q<180, original era/category selections; CR includes105-135 and VR',
            weighting='Full nominal weight including DY ParticleNet; MC sumw2; no alpha or reweighting',
            luminosity_display=138,luminosity_normalization='Original era cards, total137.64fb^-1; display rounded to138 by request',
            SR_Data_read=False,source_files=audit,validation_report_sha256=hashes),indent=2))
    render(archive,args.output)
    print('Completed 16 individual plots and 4 overview pages; SR Data never opened',flush=True)


if __name__=='__main__': main()
