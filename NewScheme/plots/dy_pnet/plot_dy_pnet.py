"""Paired-event DY ParticleNet comparisons, with correlated statistical errors."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
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
import uproot

ERAS=('2016preAPV','2016postAPV','2017','2018')
TAGS=('untag','btag','vbftag')
LABELS={'untag':'untagged','btag':'b-tag','vbftag':'VBF-tag'}
LUMI={'2016preAPV':19.52,'2016postAPV':16.81,'2017':41.48,'2018':59.83}
MASS=np.array([300,400,500,600,700,800,900,1000,1200,1400,1600,2000,2500,3000,np.inf])


def histograms(source,era,category,tag,kind):
    edges=np.arange(40.,185.,5.) if kind=='mqq' else MASS
    h=np.zeros((5,len(edges)-1)); n=0
    regions=('CR','SR') if kind=='mqq' else (kind,)
    for channel in ('ee','mumu'):
        for region in regions:
            path=source/f'DY_{channel}_{category}_{tag}_{era}_{region}.root'
            with uproot.open(path) as handle:
                for a in handle['tree'].iterate(['mass2q','m2l2q','weight','weight_no_pnet','pnet_sf'],step_size=100000,library='np'):
                    np.testing.assert_allclose(a['weight'],a['weight_no_pnet']*a['pnet_sf'],rtol=1e-14)
                    if category=='resolved': np.testing.assert_array_equal(a['weight'],a['weight_no_pnet'])
                    mask=(a['m2l2q']>=300)&(a['mass2q']>=40)&(a['mass2q']<180)
                    n+=int(mask.sum())
                    x=a['mass2q'][mask] if kind=='mqq' else a['m2l2q'][mask]
                    w=a['weight'][mask]; u=a['weight_no_pnet'][mask]
                    for i,values in enumerate((u,w,u*u,w*w,u*w)):
                        h[i]+=np.histogram(x,bins=edges,weights=values)[0]
    return edges,h,n


def ratio(h):
    u,w,vu,vw,cov=h
    valid=u>0
    r=np.full_like(u,np.nan); v=np.full_like(u,np.nan)
    r[valid]=w[valid]/u[valid]
    v[valid]=(vw[valid]+r[valid]**2*vu[valid]-2*r[valid]*cov[valid])/u[valid]**2
    np.testing.assert_array_less(-1e-10*np.maximum(1.,vw[valid]/u[valid]**2),v[valid]+1e-20)
    return r,np.sqrt(np.maximum(v,0))


def draw(edges,h,era,category,tag,kind,axes=None):
    if axes is None:
        fig,(ax,lower)=plt.subplots(2,1,figsize=(11,10),sharex=True,gridspec_kw={'height_ratios':[3.3,1],'hspace':.05})
        fig.subplots_adjust(left=.18,right=.97,top=.90,bottom=.16)
    else:
        ax,lower=axes; fig=ax.figure
    display=edges.copy()
    if np.isinf(display[-1]): display[-1]=4000
    widths=np.diff(display); scale=np.ones_like(widths) if kind=='mqq' else 100/widths
    for i,(label,color) in enumerate((('DY without ParticleNet SF','#2765a3'),('DY with ParticleNet SF','#d1493f'))):
        y=h[i]*scale; err=np.sqrt(h[i+2])*scale
        ax.stairs(y,display,color=color,lw=2.5,label=label)
        ax.fill_between(display,np.r_[y-err,(y-err)[-1]],np.r_[y+err,(y+err)[-1]],step='post',color=color,alpha=.13,lw=0)
    density=h[:2]*scale
    positive=density[density>0]
    ax.set_yscale('log'); ax.set_ylim(max(.001,float(positive.min())*.2) if len(positive) else .001,max(1.,float((h[:2]*scale).max()))*80)
    hep.cms.label('Preliminary',data=False,lumi=LUMI[era],com=13,ax=ax,fontsize=27)
    ax.text(.04,.95,f'{era}, {category.capitalize()}\n{LABELS[tag]}, '+r'$ee+\mu\mu$'+f'\n'+('CR + SR' if kind=='mqq' else kind),transform=ax.transAxes,va='top',fontsize=23)
    ax.legend(loc='upper right',fontsize=21,frameon=False,handlelength=1.5)
    ax.set_ylabel('Events / 5 GeV' if kind=='mqq' else 'Events / 100 GeV',fontsize=32)
    r,e=ratio(h); centers=(display[:-1]+display[1:])/2; valid=np.isfinite(r)
    lower.axhline(1,color='.4',lw=1)
    lower.errorbar(centers[valid],r[valid],yerr=e[valid],fmt='o',color='#d1493f',ms=5)
    vals=np.r_[r[valid]-e[valid],r[valid]+e[valid],1.]
    lo,hi=float(vals.min()),float(vals.max()); pad=max(.08,(hi-lo)*.15)
    lower.set_ylim(min(.9,lo-pad),max(1.1,hi+pad))
    lower.set_ylabel('With / without',fontsize=25)
    lower.set_xlabel((r'$m_{jj}$ [GeV]' if category=='resolved' else r'$m_J$ [GeV]') if kind=='mqq' else r'$m_{2\ell2q}$ [GeV]',fontsize=34)
    lower.set_xlim(display[0],display[-1]); lower.set_xticks([40,70,105,140,180] if kind=='mqq' else [300,1000,2000,3000,4000])
    if kind!='mqq': lower.set_xticklabels(['300','1000','2000','3000','4000']); ax.text(.98,.48,r'Last bin: $\geq 3000$ GeV',transform=ax.transAxes,ha='right',fontsize=17)
    for axis in (ax,lower): axis.tick_params(direction='in',top=True,right=True,labelsize=23); axis.grid(False)
    ax.tick_params(labelbottom=False)
    return fig


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True); parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--histograms',type=Path,help='Replot a saved archive without rereading ROOT trees')
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    plt.style.use(hep.style.CMS)
    cache={}; rows=[]; archive={}
    saved=np.load(args.histograms) if args.histograms else None
    if saved is not None:
        with args.histograms.with_name('yields.csv').open() as handle:
            saved_entries={f"{r['era']}_{r['category']}_{r['tag']}_{r['distribution']}":int(r['entries']) for r in csv.DictReader(handle)}
    for era in ERAS:
        for category in ('resolved','merged'):
            report=json.loads((args.input/f'validation_{era}_{category}.json').read_text())
            assert report['status']=='passed' and not report['pilot']
            assert (args.input/f'{era}_{category}.done').is_file()
            for tag in TAGS:
                for kind in ('mqq','CR','SR','VR'):
                    key=f'{era}_{category}_{tag}_{kind}'
                    if saved is None: edges,h,n=histograms(args.input,era,category,tag,kind)
                    else: edges,h,n=saved[key+'__edges'],saved[key+'__hist'],saved_entries[key]
                    cache[key]=(edges,h); archive[key+'__edges']=edges; archive[key+'__hist']=h
                    u,w,vu,vw,cov=h.sum(axis=1); r=w/u if u>0 else np.nan
                    err=np.sqrt(max(0.,vw+r*r*vu-2*r*cov))/u if u>0 else np.nan
                    rows.append(dict(era=era,category=category,tag=tag,distribution=kind,entries=n,
                                     without_sf=u,with_sf=w,ratio=r,ratio_stat_error=err))
        with PdfPages(args.output/f'all_{era}_DY_ParticleNet.pdf') as bundle:
            for kind in ('mqq','CR','SR','VR'):
                fig=plt.figure(figsize=(33,20)); outer=fig.add_gridspec(2,3,left=.065,right=.99,bottom=.09,top=.95,hspace=.40,wspace=.35)
                for row,category in enumerate(('resolved','merged')):
                    for col,tag in enumerate(TAGS):
                        key=f'{era}_{category}_{tag}_{kind}'; edges,h=cache[key]
                        individual=draw(edges,h,era,category,tag,kind)
                        for ext in ('png','pdf'): individual.savefig(args.output/f'{key}.{ext}',dpi=140)
                        bundle.savefig(individual); plt.close(individual)
                        sub=outer[row,col].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05)
                        ax=fig.add_subplot(sub[0]); lower=fig.add_subplot(sub[1],sharex=ax)
                        draw(edges,h,era,category,tag,kind,(ax,lower))
                for ext in ('png','pdf'): fig.savefig(args.output/f'six_categories_{era}_{kind}_DY_ParticleNet.{ext}',dpi=120)
                plt.close(fig)
        print(f'Completed plots for {era}',flush=True)
    if saved is not None: saved.close()
    np.savez_compressed(args.output/'histograms.npz',**archive)
    with (args.output/'yields.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    (args.output/'metadata.json').write_text(json.dumps(dict(input=str(args.input),selection='m2l2q>=300, 40<=m2q<180; combined ee+mumu; CR+SR union excludes VR duplication',
        uncertainties='MC statistical only, sumw2 and paired-event covariance sum(weight*weight_no_pnet); no SF systematics',
        mass_tail='All masses >=3000 folded into final displayed 3000-4000 bin',histogram_rows=['without','with','variance_without','variance_with','covariance']),indent=2))


if __name__=='__main__': main()
