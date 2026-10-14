"""User-style signal overlays, preserving old filenames in new folders."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import argparse,hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import mplhep as hep
import numpy as np
import plot_run2_kd as kd
import plot_run2_pnet as scores
import plot_fullrange_pnet as full
from prepare_histograms import MC,TAGS,CHANNELS,SIGNALS,VARS

COLORS={500:'#ff0000',1000:'#a00000',3000:'#ff88ff'}
BG_COLORS=kd.COLORS; BG_LABELS=kd.LABELS

def read_npz(path):
    with np.load(path) as f: return {k:f[k].copy() for k in f.files}

def legend(ax,background_labels):
    handles,labels=ax.get_legend_handles_labels(); lookup=dict(zip(labels,handles))
    left=[f'gg X({m})' for m in (500,1000,3000)]+background_labels[::2]
    right=[f'VBF X({m})' for m in (500,1000,3000)]+background_labels[1::2]
    n=max(len(left),len(right)); left+=['']*(n-len(left)); right+=['']*(n-len(right))
    labels=left+right; blank=Line2D([],[],color='none')
    ax.legend([lookup.get(k,blank) for k in labels],labels,ncol=2,loc='upper right',fontsize=20,
              frameon=False,handlelength=1.5,columnspacing=.9,labelspacing=.28,
              title=r'Signals: $M_X$ [GeV], $\sigma=1$ pb each',title_fontsize=17)

def signals(ax,archive,key,edges,scale=1.):
    maximum=0.
    for sample in SIGNALS:
        mode=sample[:3]; mass=int(sample[3:]); values=archive[key+'__'+sample][0]*scale
        label=('gg' if mode=='ggh' else 'VBF')+f' X({mass})'
        ax.stairs(values,edges,color=COLORS[mass],ls=':' if mode=='ggh' else '-',lw=2.3,label=label,zorder=4)
        maximum=max(maximum,float(values.max(initial=0)))
    return maximum

def draw_overlay(module,bg,sig,cat,region,var,linear=False,axes=None):
    fig=module.draw(bg,cat,region,var,axes)
    ax=axes[0] if axes is not None else fig.axes[0]
    old_labels=[t.get_text() for t in ax.get_legend().get_texts()]; ax.get_legend().remove()
    key=f'{cat}_{region}_{var}'; edges=bg[key+'__edges']
    maximum=signals(ax,sig,key,edges)
    total=sum(bg[key+'__'+p][0] for p in MC); error=np.sqrt(sum(bg[key+'__'+p][1] for p in MC))
    maximum=max(maximum,float((total+error).max()))
    if key+'__DY_no_pnet' in bg:
        off=total-bg[key+'__DY'][0]+bg[key+'__DY_no_pnet'][0]
        maximum=max(maximum,float(off.max()))
    if region=='CR':
        d=bg[key+'__Data'][0]; maximum=max(maximum,float((d+module.poisson(d)[1]).max()))
    if linear: ax.set_yscale('linear'); ax.set_ylim(0,max(maximum,1.)*2.6)
    else:
        candidates=[total]+[sig[key+'__'+p][0] for p in SIGNALS]
        positive=np.concatenate([v[v>0] for v in candidates])
        low=max(.001,float(positive.min())*.2) if len(positive) else .001
        ax.set_yscale('log'); ax.set_ylim(low,max(maximum,1.)*300)
    legend(ax,old_labels)
    return fig

def save(fig,out,name,bundle=None):
    for ext in ('png','pdf'): fig.savefig(out/f'{name}.{ext}',dpi=140 if ext=='png' else None)
    if bundle is not None: bundle.savefig(fig)

def render_kd(bg,sig,out):
    out.mkdir(parents=True,exist_ok=True)
    with PdfPages(out/'all_Run2_discriminants.pdf') as bundle:
        for region in ('SR','CR'):
            for cat in ('resolved','merged'):
                overview=plt.figure(figsize=(33,10)); grid=overview.add_gridspec(1,3,left=.06,right=.985,bottom=.16,top=.90,wspace=.32)
                for i,var in enumerate(('kd','kd_jjvbf','kd_jvbf')):
                    fig=draw_overlay(kd,bg,sig,cat,region,var); save(fig,out,f'{cat}_{region}_{var}',bundle); plt.close(fig)
                    if region=='CR':
                        sub=grid[i].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05); ax=overview.add_subplot(sub[0]); ratio=overview.add_subplot(sub[1],sharex=ax)
                    else: ax=overview.add_subplot(grid[i]); ratio=None
                    draw_overlay(kd,bg,sig,cat,region,var,axes=(ax,ratio))
                save(overview,out,f'{cat}_{region}_three_discriminants'); plt.close(overview)

def render_pnet(bg,sig,out,linear=False,fullrange=False):
    out.mkdir(parents=True,exist_ok=True)
    module=full if fullrange else scores
    bundle_name='all_Run2_fullrange_PNet.pdf' if fullrange else 'all_Run2_PNet_scores.pdf'
    with PdfPages(out/bundle_name) as bundle:
        if fullrange:
            combinations=[('SR','particleNetZvsQCD'),('CR','particleNetZvsQCD')]
            overview=plt.figure(figsize=(22,10)); grid=overview.add_gridspec(1,2,left=.095,right=.985,bottom=.16,top=.90,wspace=.32)
            for i,(region,var) in enumerate(combinations):
                fig=draw_overlay(module,bg,sig,'merged',region,var,linear); save(fig,out,f'merged_{region}_particleNetZvsQCD_fullrange',bundle); plt.close(fig)
                if region=='CR':
                    sub=grid[i].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05); ax=overview.add_subplot(sub[0]); ratio=overview.add_subplot(sub[1],sharex=ax)
                else: ax=overview.add_subplot(grid[i]); ratio=None
                draw_overlay(module,bg,sig,'merged',region,var,linear,axes=(ax,ratio))
            save(overview,out,'SR_CR_fullrange_particleNetZvsQCD'); plt.close(overview)
        else:
            for region in ('SR','CR'):
                overview=plt.figure(figsize=(22,10)); grid=overview.add_gridspec(1,2,left=.095,right=.985,bottom=.16,top=.90,wspace=.32)
                for i,var in enumerate(('particleNetZvsQCD','particleNetZbbvslight')):
                    fig=draw_overlay(module,bg,sig,'merged',region,var,linear); save(fig,out,f'merged_{region}_{var}',bundle); plt.close(fig)
                    if region=='CR':
                        sub=grid[i].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05); ax=overview.add_subplot(sub[0]); ratio=overview.add_subplot(sub[1],sharex=ax)
                    else: ax=overview.add_subplot(grid[i]); ratio=None
                    draw_overlay(module,bg,sig,'merged',region,var,linear,axes=(ax,ratio))
                save(overview,out,f'merged_{region}_PNet_scores'); plt.close(overview)

def mass_parts(archive,cat,tag,ch):
    tags=TAGS if tag=='all' else (tag,); channels=CHANNELS if ch=='combined' else (ch,)
    return {p:sum((archive[f'{cat}__{t}__{c}__{p}'] for t in tags for c in channels),np.zeros((2,len(archive['mass_edges'])-1))) for p in MC+SIGNALS}

def draw_mass(archive,cat,tag,ch,ax=None):
    if ax is None:
        fig,ax=plt.subplots(figsize=(11,9)); fig.subplots_adjust(left=.18,right=.97,top=.90,bottom=.16)
    else: fig=ax.figure
    edges=archive['mass_edges'].copy(); edges[-1]=4000.; scale=50./np.diff(edges)
    parts=mass_parts(archive,cat,tag,ch); total=np.zeros(len(scale)); variance=total.copy()
    for p in MC:
        values=parts[p][0]*scale
        ax.bar(edges[:-1],values,width=np.diff(edges),bottom=total,align='edge',color=BG_COLORS[p],edgecolor='none',label=BG_LABELS[p]); total+=values; variance+=parts[p][1]*scale**2
    err=np.sqrt(variance); step=lambda a:np.r_[a,a[-1]]; band=r'Bkg. MC stat. ($1\sigma$)'
    ax.fill_between(edges,step(total-err),step(total+err),step='post',color='#4c566a',alpha=.3,lw=0,label=band)
    signal_archive={'mass__'+p:parts[p] for p in SIGNALS}
    maximum=max(float((total+err).max()),signals(ax,signal_archive,'mass',edges,scale))
    ax.set_yscale('log'); ax.set_ylim(.001,max(maximum,1.)*300)
    hep.cms.label('Preliminary',data=False,lumi=138,com=13,ax=ax,fontsize=27)
    taglabel={'all':'All tags','untag':'Untagged','btag':'b-tag','vbftag':'VBF-tag'}[tag]
    chlabel={'ee':r'$ee$','mumu':r'$\mu\mu$','combined':r'$ee+\mu\mu$'}[ch]
    ax.text(.035,.95,f'Run 2, {cat.capitalize()}\n{taglabel}\n'+chlabel+'\nSR (blinded)',transform=ax.transAxes,va='top',fontsize=21)
    legend(ax,[BG_LABELS[p] for p in MC]+[band])
    ax.set_ylabel('Events / 50 GeV',fontsize=32); ax.set_xlabel(r'$m_{\ell\ell jj}$ [GeV]' if cat=='resolved' else r'$m_{\ell\ell J}$ [GeV]',fontsize=34)
    ax.set_xlim(300,4000); ax.set_xticks([500,1000,2000,3000,3750]); ax.set_xticklabels(['500','1000','2000','3000',r'$\geq3500$'])
    ax.tick_params(direction='in',top=True,right=True,labelsize=23); ax.grid(False)
    ax.get_xticklabels()[-1].set_horizontalalignment('right')
    return fig

def render_mass(archive,out):
    out.mkdir(parents=True,exist_ok=True)
    with PdfPages(out/'all_SR_m2l2q_signals.pdf') as bundle:
        for cat in ('resolved','merged'):
            for tag in TAGS+('all',):
                for ch in CHANNELS+('combined',):
                    fig=draw_mass(archive,cat,tag,ch); save(fig,out,f'{cat}_{tag}_{ch}_SR_m2l2q',bundle); plt.close(fig)
        for ch in CHANNELS+('combined',):
            fig=plt.figure(figsize=(33,20)); grid=fig.add_gridspec(2,3,left=.07,right=.985,bottom=.09,top=.95,wspace=.34,hspace=.36)
            for i,cat in enumerate(('resolved','merged')):
                for j,tag in enumerate(TAGS): draw_mass(archive,cat,tag,ch,fig.add_subplot(grid[i,j]))
            save(fig,out,f'six_categories_SR_{ch}_m2l2q'); plt.close(fig)
        fig=plt.figure(figsize=(22,10)); grid=fig.add_gridspec(1,2,left=.095,right=.985,bottom=.16,top=.90,wspace=.32)
        for i,cat in enumerate(('resolved','merged')): draw_mass(archive,cat,'all','combined',fig.add_subplot(grid[i]))
        save(fig,out,'inclusive_SR_m2l2q'); plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--analysis',type=Path,required=True); p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--mass-only',action='store_true',help='Rerender only the new SR mass plots from cache.')
    a=p.parse_args(); plot=a.analysis/'Plot'; plt.style.use(hep.style.CMS)
    if a.mass_only:
        render_mass(read_npz(a.cache/'mass_histograms.npz'),plot/'run2_SR_m2l2q_signals')
        return
    sig=read_npz(a.cache/'signal_histograms.npz'); rawsig=read_npz(a.cache/'signal_fullrange_histograms.npz'); mass=read_npz(a.cache/'mass_histograms.npz')
    kd_bg=read_npz(plot/'run2_discriminants_PNet'/'histograms.npz'); score_bg=read_npz(plot/'run2_PNet_scores'/'histograms.npz'); full_bg=read_npz(plot/'run2_PNet_fullrange'/'histograms.npz')
    render_kd(kd_bg,sig,plot/'run2_discriminants_PNet_signals'); print('Rendered MELA signal overlays',flush=True)
    for linear in (False,True):
        render_pnet(score_bg,sig,plot/('run2_PNet_scores'+('_linear' if linear else '')+'_signals'),linear)
        render_pnet(full_bg,rawsig,plot/('run2_PNet_fullrange'+('_linear' if linear else '')+'_signals'),linear,True)
        print('Rendered ParticleNet signal overlays',linear,flush=True)
    render_mass(mass,plot/'run2_SR_m2l2q_signals'); print('Rendered all blinded SR mass plots',flush=True)
    names=['run2_discriminants_PNet_signals','run2_PNet_scores_signals','run2_PNet_scores_linear_signals','run2_PNet_fullrange_signals','run2_PNet_fullrange_linear_signals','run2_SR_m2l2q_signals']
    for name in names:
        out=plot/name
        (out/'README.md').write_text('# Run-2 signal overlays\n\n'+
            'Separate signal-overlay version; original plot folders unchanged. Backgrounds/CR data and statistical bands are unchanged cached values. SR Data is not read or plotted. Signals are six independent 1 pb benchmark hypotheses, not stacked or added to ratio denominators.\n\n'+
            'Style: ggF dotted, VBF solid; 500 red, 1000 dark red, 3000 pink, with paired legend columns. Nominal corrected signal weights for MELA, selected ParticleNet and mass plots; no ParticleNet SF for full-range ZvsQCD. Selected-score plots retain the DY-SF-off total overlay, switching only DY.\n\n'+
            'All eras, original luminosity normalization (138 fb^-1 display). Direct DY MC for new SR mass plots. Mass bins share an event-filled grid across channels/categories:50 GeV300-1500,100 GeV1500-3500,>=3500 retained in a displayed3500-4000 overflow bin and density-scaled to Events/50GeV. No signal-free high-mass CR is assumed.\n\n'+
            'Scripts, histogram caches, yields and source provenance: ../signal_overlay_Run2_20261013/. Individual PNG/PDF files and overview pages are provided.\n')
    (a.cache/'outputs.json').write_text(json.dumps({name:dict(png=len(list((plot/name).glob('*.png'))),pdf=len(list((plot/name).glob('*.pdf')))) for name in names},indent=2))
if __name__=='__main__': main()
