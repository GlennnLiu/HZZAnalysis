"""Measured in-frame legend clearance and compact limits; CR mass plots."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='1'
import argparse,csv,json,hashlib,importlib.util
from pathlib import Path
import numpy as np
import uproot
import render_plots as r
from prepare_histograms import ERAS,TAGS,CHANNELS,MC,MASS_EDGES,fill
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import StepPatch,Rectangle

AUDIT=[]
def layout(fig):
    if getattr(fig,'_revised',False): return
    fig._revised=True
    for ax in fig.axes:
        leg=ax.get_legend()
        if leg is None: continue
        labels=[t.get_text() for t in leg.get_texts()]
        title=leg.get_title().get_text()
        handles,artistlabels=ax.get_legend_handles_labels(); lookup=dict(zip(artistlabels,handles))
        labels=[s.replace('DY + jets, SF on','DY + jets').replace('Total, DY SF off','DY SF off (total)') for s in labels]
        lookup={k.replace('DY + jets, SF on','DY + jets').replace('Total, DY SF off','DY SF off (total)'):v for k,v in lookup.items()}
        signal=any('X(' in s for s in labels)
        region=None
        for t in ax.texts:
            if t.get_text().startswith('Run 2'):
                region=t; t.set_text(t.get_text().replace('All tags, ','All tags\n')); t.set_fontsize(20)
            elif 'Last bin:' in t.get_text():
                t.set_visible(False)
        leg.remove()
        spans=[]
        for p in ax.patches:
            if isinstance(p,StepPatch):
                data=p.get_data(); spans.extend(zip(data.edges[:-1],data.edges[1:],data.values))
            elif isinstance(p,Rectangle): spans.append((p.get_x(),p.get_x()+p.get_width(),p.get_y()+p.get_height()))
        for line in ax.lines:
            spans.extend((x,x,y) for x,y in zip(line.get_xdata(),line.get_ydata()))
        for collection in ax.collections:
            for path in collection.get_paths(): spans.extend((x,x,y) for x,y in path.vertices)
        spans=np.asarray(spans,dtype=float); spans=spans[np.isfinite(spans).all(axis=1)&(spans[:,2]>0)]
        positives=spans[:,2]
        ymax=max(float(ax.dataLim.ymax),float(positives.max()) if len(positives) else 1.)
        logarithmic=ax.get_yscale()=='log'; low=max(.001,float(positives.min())*.5) if len(positives) else .001
        choices=[]
        for side in ('right','left'):
            if ax.get_legend() is not None: ax.get_legend().remove()
            new=ax.legend([lookup.get(s,Line2D([],[],color='none')) for s in labels],labels,
                          ncol=2 if signal else 1,fontsize=19 if signal else 20,frameon=False,
                          handlelength=1.2,columnspacing=.7,labelspacing=.23,
                          loc='upper '+side,bbox_to_anchor=(.98 if side=='right' else .02,.975),
                          title=title,title_fontsize=16,borderaxespad=0)
            if region is not None:
                region.set_position((.025 if side=='right' else .975,.97)); region.set_ha('left' if side=='right' else 'right')
            ax.set_ylim(low,max(ymax*1.25,low*10)) if logarithmic else ax.set_ylim(0,max(ymax,1.)*1.15)
            fig.canvas.draw(); renderer=fig.canvas.get_renderer()
            boxes=[new.get_window_extent(renderer)]+([region.get_window_extent(renderer)] if region else [])
            overlap=len(boxes)>1 and boxes[0].overlaps(boxes[1])
            required=np.log(max(ymax*1.25,low*10)) if logarithmic else max(ymax,1.)*1.15
            for box in boxes:
                xy=ax.transAxes.inverted().transform(box.get_points()); fraction=max(.12,xy[0,1]-.025)
                x1,x2=ax.transData.inverted().transform([[box.x0,ax.bbox.y0],[box.x1,ax.bbox.y0]])[:,0]
                masked=spans[(spans[:,0]<=x2)&(spans[:,1]>=x1)]
                if len(masked):
                    v=masked[:,2].max()
                    required=max(required,(np.log(v)-(1-fraction)*np.log(low))/fraction if logarithmic else v/fraction)
            choices.append((required if not overlap else float('inf'),side))
        required,side=min(choices)
        assert np.isfinite(required),'Legend and category text collide horizontally'
        if side=='right':
            new.set_bbox_to_anchor((.98,.975)); new._loc=1
            if region: region.set_position((.025,.97)); region.set_ha('left')
        ax.set_ylim(low,np.exp(required)*1.02) if logarithmic else ax.set_ylim(0,required*1.02)
        AUDIT.append(dict(labels=labels,ylim=list(ax.get_ylim()),yscale=ax.get_yscale(),legend_inside_plot=True,legend_side=side,measured_clearance=True))

original_save=r.save
def save(fig,*args,**kwargs): layout(fig); return original_save(fig,*args,**kwargs)
r.save=save

def load_cr(analysis,cache):
    path=cache/'CR_mass_histograms.npz'
    if path.exists(): return r.read_npz(path)
    root=analysis/'FitTemplate'/'RedTrees_20261009_pnet'; store={}; rows=[]
    for cat in ('resolved','merged'):
        for era in ERAS:
            for tag in TAGS:
                for ch in CHANNELS:
                    for proc in MC+('Data',):
                        source=root/f'{proc}_{ch}_{cat}_{tag}_{era}_CR.root'
                        sums=np.zeros(2); count=0
                        with uproot.open(source) as f:
                            for v in f['tree'].iterate(['m2l2q','mass2q','weight'],step_size=100000,library='np'):
                                q=v['mass2q']; selected=(v['m2l2q']>=300)&(q>=40)&(q<180)
                                assert not ((q[selected]>70)&(q[selected]<105)).any()
                                w=v['weight'][selected].astype(float)
                                assert np.isfinite(w).all()
                                if proc=='Data': np.testing.assert_array_equal(w,np.ones(len(w)))
                                fill(store,f'{cat}__{tag}__{ch}__{proc}',v['m2l2q'][selected],w,MASS_EDGES)
                                sums+=[w.sum(),(w*w).sum()]; count+=len(w)
                        rows.append(dict(file=str(source),entries=count,sumw=sums[0],sumw2=sums[1]))
        print('Loaded CR mass:',cat,flush=True)
    for cat in ('resolved','merged'):
        for tag in TAGS:
            for ch in CHANNELS:
                for proc in MC+('Data',):
                    suffix=f'{proc}_{ch}_{cat}_{tag}_'
                    wanted=[a for a in rows if Path(a['file']).name.startswith(suffix)]
                    np.testing.assert_allclose(store[f'{cat}__{tag}__{ch}__{proc}'].sum(axis=1),[sum(a[k] for a in wanted) for k in ('sumw','sumw2')],rtol=1e-10,atol=1e-6)
    np.savez_compressed(path,mass_edges=MASS_EDGES,**store)
    with (cache/'CR_mass_input_yields.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    return dict(mass_edges=MASS_EDGES,**store)

def draw_cr(archive,cat,tag,ch,axes=None):
    if axes is None:
        fig,(ax,ratio)=r.plt.subplots(2,1,figsize=(11,10),sharex=True,gridspec_kw={'height_ratios':[3.3,1],'hspace':.05})
        fig.subplots_adjust(left=.18,right=.97,top=.90,bottom=.16)
    else: ax,ratio=axes; fig=ax.figure
    tags=TAGS if tag=='all' else (tag,); channels=CHANNELS if ch=='combined' else (ch,)
    parts={p:sum((archive[f'{cat}__{t}__{c}__{p}'] for t in tags for c in channels)) for p in MC+('Data',)}
    edges=archive['mass_edges'].copy(); edges[-1]=4000.; scale=50/np.diff(edges); centers=(edges[:-1]+edges[1:])/2
    total=np.zeros(len(scale)); variance=total.copy()
    for p in MC:
        h=parts[p]; ax.bar(edges[:-1],h[0]*scale,width=np.diff(edges),bottom=total*scale,align='edge',color=r.BG_COLORS[p],edgecolor='none',label=r.BG_LABELS[p]); total+=h[0]; variance+=h[1]
    error=np.sqrt(variance); step=lambda v:np.r_[v,v[-1]]
    ax.fill_between(edges,step((total-error)*scale),step((total+error)*scale),step='post',color='#4c566a',alpha=.3,lw=0,label=r'Bkg. MC stat. ($1\sigma$)')
    data=parts['Data'][0]; de=r.kd.poisson(data)
    ax.errorbar(centers,data*scale,yerr=de*scale,fmt='o',color='black',ms=5,lw=1.2,label='Data',zorder=5)
    good=total>0; rel=np.divide(error,total,out=np.full_like(total,np.nan),where=good)
    ratio.fill_between(edges,step(1-rel),step(1+rel),step='post',color='#4c566a',alpha=.3,lw=0)
    ratio.errorbar(centers[good],data[good]/total[good],yerr=de[:,good]/total[good],fmt='o',color='black',ms=5,lw=1.2)
    ratio.axhline(1,color='.4',lw=1); ratio.set_ylim(0,2.5); ratio.set_ylabel('Data / MC',fontsize=27)
    beyond=good&(data>2.5*total); ratio.plot(centers[beyond],np.full(beyond.sum(),2.43),'^',color='black',ms=7)
    ax.set_yscale('log'); r.hep.cms.label('Preliminary',data=True,lumi=138,com=13,ax=ax,fontsize=27)
    taglabel={'all':'All tags','untag':'Untagged','btag':'b-tag','vbftag':'VBF-tag'}[tag]
    chlabel={'ee':r'$ee$','mumu':r'$\mu\mu$','combined':r'$ee+\mu\mu$'}[ch]
    ax.text(.035,.95,f'Run 2, {cat.capitalize()}\n{taglabel}\n'+chlabel+'\nCR',transform=ax.transAxes,va='top',fontsize=21)
    ax.legend(loc='upper right',fontsize=20,frameon=False)
    ax.set_ylabel('Events / 50 GeV',fontsize=32); ratio.set_xlabel(r'$m_{\ell\ell jj}$ [GeV]' if cat=='resolved' else r'$m_{\ell\ell J}$ [GeV]',fontsize=34)
    ratio.set_xlim(300,4000); ratio.set_xticks([500,1000,2000,3000,3750]); ratio.set_xticklabels(['500','1000','2000','3000',r'$\geq3500$']); ratio.get_xticklabels()[-1].set_horizontalalignment('right')
    for a in (ax,ratio): a.tick_params(direction='in',top=True,right=True,labelsize=23); a.grid(False)
    ax.tick_params(labelbottom=False)
    return fig

def render_cr(archive,out):
    out.mkdir(exist_ok=True)
    with PdfPages(out/'all_CR_m2l2q_DataMC.pdf') as bundle:
        for cat in ('resolved','merged'):
            for tag in TAGS+('all',):
                for ch in CHANNELS+('combined',):
                    fig=draw_cr(archive,cat,tag,ch); save(fig,out,f'{cat}_{tag}_{ch}_CR_m2l2q',bundle); r.plt.close(fig)
        for ch in CHANNELS+('combined',):
            fig=r.plt.figure(figsize=(33,22)); grid=fig.add_gridspec(2,3,left=.07,right=.985,bottom=.09,top=.95,wspace=.34,hspace=.40)
            for i,cat in enumerate(('resolved','merged')):
                for j,tag in enumerate(TAGS):
                    sub=grid[i,j].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05); ax=fig.add_subplot(sub[0]); ratio=fig.add_subplot(sub[1],sharex=ax)
                    draw_cr(archive,cat,tag,ch,(ax,ratio))
            save(fig,out,f'six_categories_CR_{ch}_m2l2q'); r.plt.close(fig)
        fig=r.plt.figure(figsize=(22,11)); grid=fig.add_gridspec(1,2,left=.095,right=.985,bottom=.16,top=.90,wspace=.32)
        for i,cat in enumerate(('resolved','merged')):
            sub=grid[i].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05); ax=fig.add_subplot(sub[0]); ratio=fig.add_subplot(sub[1],sharex=ax); draw_cr(archive,cat,'all','combined',(ax,ratio))
        save(fig,out,'inclusive_CR_m2l2q'); r.plt.close(fig)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--analysis',type=Path,required=True); p.add_argument('--preview',action='store_true'); p.add_argument('--part',choices=('all','signals','backgrounds','cr'),default='all'); a=p.parse_args()
    plot=a.analysis/'Plot'; cache=plot/'signal_overlay_Run2_20261013'; r.plt.style.use(r.hep.style.CMS)
    bg=r.read_npz(plot/'run2_discriminants_PNet/histograms.npz'); scores=r.read_npz(plot/'run2_PNet_scores/histograms.npz'); full=r.read_npz(plot/'run2_PNet_fullrange/histograms.npz')
    sig=r.read_npz(cache/'signal_histograms.npz'); rawsig=r.read_npz(cache/'signal_fullrange_histograms.npz'); mass=r.read_npz(cache/'mass_histograms.npz')
    if a.preview:
        out=cache/'layout_preview'; out.mkdir(exist_ok=True)
        fig=r.draw_overlay(r.scores,scores,sig,'merged','SR','particleNetZvsQCD'); save(fig,out,'pnet_signal'); r.plt.close(fig)
        fig=r.draw_mass(mass,'merged','vbftag','combined'); save(fig,out,'mass_signal'); r.plt.close(fig)
        fig=r.scores.draw(scores,'merged','CR','particleNetZvsQCD'); save(fig,out,'pnet_background'); r.plt.close(fig)
        return
    if a.part in ('all','cr'):
        cr=load_cr(a.analysis,cache); render_cr(cr,plot/'run2_CR_m2l2q_DataMC'); print('Completed CR mass plots, no signals',flush=True)
    if a.part in ('all','signals'):
        overlay=r.draw_overlay
        def sr_signals_only(mod,b,s,cat,reg,var,linear=False,axes=None):
            if reg=='SR': return overlay(mod,b,s,cat,reg,var,linear,axes)
            fig=mod.draw(b,cat,reg,var,axes)
            if linear: (axes[0] if axes is not None else fig.axes[0]).set_yscale('linear')
            return fig
        r.draw_overlay=sr_signals_only
        r.render_kd(bg,sig,plot/'run2_discriminants_PNet_signals_layout'); print('Revised MELA signals',flush=True)
        for linear in (False,True):
            suffix=('_linear' if linear else '')+'_signals_layout'
            r.render_pnet(scores,sig,plot/('run2_PNet_scores'+suffix),linear)
            r.render_pnet(full,rawsig,plot/('run2_PNet_fullrange'+suffix),linear,True)
        r.render_mass(mass,plot/'run2_SR_m2l2q_signals_layout'); print('Revised signal plots',flush=True)
    if a.part in ('signals','cr'):
        (cache/f'layout_validation_{a.part}.json').write_text(json.dumps(AUDIT,indent=2)); return
    # Use exactly the same background caches and renderers, applying only layout.
    for module,name,archive,method in ((r.kd,'run2_discriminants_PNet',bg,'kd'),(r.scores,'run2_PNet_scores',scores,'scores'),(r.full,'run2_PNet_fullrange',full,'full')):
        original=module.draw
        # r rendering invokes draw per axis. Layout is applied once at save, not here.
        r.draw_overlay=lambda mod,b,s,cat,reg,var,linear=False,axes=None: mod.draw(b,cat,reg,var,axes)
        if method=='kd': r.render_kd(archive,sig,plot/(name+'_layout'))
        else:
            for linear in (False,True):
                if linear:
                    def linear_draw(*args,_original=original,**kwargs):
                        fig=_original(*args,**kwargs); axes=args[-1]; (axes[0] if axes is not None else fig.axes[0]).set_yscale('linear'); return fig
                    module.draw=linear_draw
                r.render_pnet(archive,sig if method=='scores' else rawsig,plot/(name+('_linear' if linear else '')+'_layout'),linear,method=='full')
                module.draw=original
        print('Revised background plots',name,flush=True)
    spec=importlib.util.spec_from_file_location('zkin',a.analysis/'NewScheme/plots/z_kinematics/plot_run2_z.py'); z=importlib.util.module_from_spec(spec); spec.loader.exec_module(z)
    original=z.draw
    # Its renderer saves directly; apply layout immediately for individual figures,
    # and when all four overview panels have been drawn.
    def zdraw(*args,**kwargs):
        fig=original(*args,**kwargs)
        if len(args)<5 or args[4] is None: layout(fig)
        elif sum(ax.get_legend() is not None for ax in fig.axes)==4: layout(fig)
        return fig
    z.draw=zdraw; zout=plot/'run2_z_kinematics_PNet_layout'; zout.mkdir(exist_ok=True); z.render(r.read_npz(plot/'run2_z_kinematics_PNet/histograms.npz'),zout)
    (cache/'layout_validation.json').write_text(json.dumps(AUDIT,indent=2))
    print('PASS: measured in-frame legend clearance; CR mass yields/variances preserved; no SR data opened.',flush=True)
if __name__=='__main__': main()
