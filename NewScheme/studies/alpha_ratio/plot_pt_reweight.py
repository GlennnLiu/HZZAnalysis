"""CMS-style PDF/PNG presentation of the isolated DY pT study."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import zipfile
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import mplhep as hep
import numpy as np

TAGS=('btag','untag','vbftag')
CHANNELS=('ee','mumu','combined')
COLORS={'nominal':'#e68a00','shared':'#2166ac','category':'#b2182b','normalization':'#666666'}
LABELS={'btag':'b-tag','untag':'untagged','vbftag':'VBF-tag'}
CH={'ee':r'$ee$','mumu':r'$\mu\mu$','combined':r'$ee+\mu\mu$'}


def display_edges(edges):
    out=edges.copy()
    if np.isinf(out[-1]):out[-1]=out[-2]+max(250,min(800,(out[-2]-out[0])*.35))
    return out


def step(v):return np.r_[v,v[-1]]


def prefix(cat,tag,ch,target,scope,cal):return '__'.join((cat,tag,ch,target,scope,cal))+'__'


def massplot(data,summary,cat,tag,ch,target,method,cal,compact=False,axes=None):
    if axes is None:
        fig,axes=plt.subplots(2,1,figsize=(10,9),sharex=True,
                            gridspec_kw={'height_ratios':[3,1],'hspace':.05})
        fig.subplots_adjust(top=.91,bottom=.16,left=.15,right=.96)
    else:fig=axes[0].figure
    ax,lower=axes
    root=prefix(cat,tag,ch,target,'shared',cal)
    edges=data[root+'edges'];x=display_edges(edges);center=(x[:-1]+x[1:])/2
    observed=data[root+'observed']
    ax.errorbar(center,observed[0],yerr=np.sqrt(observed[1]),fmt='o',color='black',
                markersize=3 if compact else 4,label='Data - TTplusWW - VZ',zorder=5)
    nominalkey=('dy_nominal' if method=='direct' else
                'alpha_inclusiveCR_nominal' if target=='SR' else 'alpha_nominal')
    reweightedkey=('dy_reweighted' if method=='direct' else
                   'alpha_inclusiveCR_reweighted' if target=='SR' else 'alpha_reweighted')
    curves=[('nominal',root+nominalkey,'Nominal DY' if method=='direct' else 'Nominal alpha'),
            ('shared',root+reweightedkey,'pT reweight: shared tags'),
            ('category',prefix(cat,tag,ch,target,'category',cal)+reweightedkey,'pT reweight: per category')]
    if method=='direct':curves.append(('normalization',prefix(cat,tag,ch,target,'category',cal)+'dy_normonly','Normalization only'))
    ymax=float(np.max(observed[0]+np.sqrt(observed[1])));ymin=min(0,float(np.min(observed[0]-np.sqrt(observed[1]))))
    for scope,key,label in curves:
        val,var=data[key];err=np.sqrt(np.maximum(var,0))
        ax.step(x,step(val),where='post',color=COLORS[scope],linewidth=1.8,
                linestyle='--' if scope in ('category','normalization') else '-',label=label)
        if scope!='normalization':
            ax.fill_between(x,step(val-err),step(val+err),step='post',color=COLORS[scope],alpha=.13,linewidth=0)
        ymax=max(ymax,float(np.max(val+err)));ymin=min(ymin,float(np.min(val-err)))
        if scope=='normalization':continue
        closure=data[key+'_closure_bins'];use=data[root+'chi2_mask'].astype(bool)&np.isfinite(closure[0])
        offset={'nominal':-.10,'shared':0,'category':.10}[scope]
        lower.errorbar(center[use]+offset*np.diff(x)[use],closure[0,use],
            yerr=np.sqrt(np.maximum(0,closure[1,use])),fmt='o',ms=2.5 if compact else 3,
            color=COLORS[scope],elinewidth=.8,capsize=0)
    ax.set_ylim(ymin-.035*(ymax-ymin),max(1,ymax)*1.6)
    hep.cms.label('Work in progress',data=True,lumi=59.8,com=13,fontsize=10 if compact else 19,ax=ax)
    ax.text(.035,.965,f'2018, {target}\n{cat.capitalize()}, {LABELS[tag]}\n{CH[ch]}',
            transform=ax.transAxes,ha='left',va='top',fontsize=10 if compact else 14)
    marker='DY' if method=='direct' else 'alpha'
    rwrow=summary[(cat,tag,ch,target,'shared',cal,reweightedkey)]
    nominalrow=summary[(cat,tag,ch,target,'shared',cal,nominalkey)]
    text=(f"Observed/predicted:\nNominal: {float(nominalrow['closure']):.3f}\npT shared: {float(rwrow['closure']):.3f}")
    ax.text(.975,.965,text,transform=ax.transAxes,ha='right',va='top',fontsize=8 if compact else 11)
    if not compact:ax.legend(loc='upper right',bbox_to_anchor=(.99,.76),fontsize=11,frameon=False)
    lower.axhline(1,color='0.3',lw=.9)
    lower.set_ylabel('Obs. / pred.',fontsize=11 if compact else 17)
    # All bins appear above and in the tables; ratio errors are drawn for stable bins only.
    visible=[]
    for scope,key,_ in curves[:3]:
        ratio=data[key+'_closure_bins'][0]
        mask=data[root+'chi2_mask'].astype(bool)&np.isfinite(ratio)
        visible.extend(ratio[mask].tolist())
    low=min(.7,min(visible)-.1) if visible else .5
    high=max(1.3,max(visible)+.1) if visible else 1.5
    lower.set_ylim(low,high)
    ticks=[v for v in (300,500,750,1100,1600,2300) if v<x[-2]]+[center[-1]]
    lower.set_xticks(ticks);lower.set_xticklabels([f'{v:g}' for v in ticks[:-1]]+[rf'$\geq {x[-2]:g}$'])
    lower.set_xlabel(r'$m_{2\ell2q}$ [GeV]' if cat=='resolved' else r'$m_{2\ell J}$ [GeV]',fontsize=13 if compact else 20)
    for axis in axes:
        axis.set_xlim(x[0],x[-1]);axis.grid(False)
        axis.tick_params(which='both',direction='in',top=True,right=True,labelsize=9 if compact else 13)
        axis.tick_params(which='major',length=5 if compact else 7)
        axis.tick_params(which='minor',length=2 if compact else 3)
    ax.tick_params(labelbottom=False)
    ax.set_ylabel('Events / bin',fontsize=12 if compact else 19)
    if not compact:
        qualifier=('VR enters pT calibration: in-sample diagnostic' if target=='VR' and cal=='inclVR'
                   else 'SR and VR held out of pT calibration' if cal=='heldoutVR' else 'SR held out of pT calibration')
        fig.text(.5,.065,qualifier,ha='center',fontsize=11)
        fig.text(.5,.025,'Statistical errors only; last bin includes overflow. Closure panel: populated bins only.',ha='center',fontsize=10)
    return fig


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    args=p.parse_args();output=args.input
    if (output/'all_pt_reweight_plots.pdf').exists():raise FileExistsError('Plots already exist')
    with np.load(output/'results.npz') as f:data={key:f[key] for key in f.files}
    with (output/'summary.csv').open() as f:rows=list(csv.DictReader(f))
    summary={tuple(row[k] for k in ('category','tag','channel','target','scope','calibration','prediction')):row for row in rows}
    plt.style.use(hep.style.CMS)
    plt.rcParams.update({'font.size':14,'axes.labelpad':7,'pdf.fonttype':42,'axes.linewidth':1.1})
    count=0
    with PdfPages(output/'all_pt_reweight_plots.pdf') as bundle:
        def save(fig,stem):
            nonlocal count
            fig.savefig(output/(stem+'.pdf'),bbox_inches='tight')
            fig.savefig(output/(stem+'.png'),dpi=180,bbox_inches='tight')
            bundle.savefig(fig,bbox_inches='tight');plt.close(fig);count+=1
        for method,target in [('direct','SR'),('alpha','SR'),('alpha','VR')]:
            fig=plt.figure(figsize=(17,11.5));outer=fig.add_gridspec(2,3,wspace=.32,hspace=.38)
            for i,(cat,tag) in enumerate((c,t) for c in ('resolved','merged') for t in TAGS):
                grid=outer[i//3,i%3].subgridspec(2,1,height_ratios=[3,1],hspace=.05)
                ax=fig.add_subplot(grid[0]);lower=fig.add_subplot(grid[1],sharex=ax)
                massplot(data,summary,cat,tag,'combined',target,method,'inclVR',compact=True,axes=(ax,lower))
            handles,labels=fig.axes[0].get_legend_handles_labels()
            fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.995),ncol=4 if method=='direct' else 3,fontsize=11,frameon=False)
            qualifier='VR used in calibration; in-sample diagnostic' if target=='VR' else 'SR held out of calibration'
            fig.text(.5,.015,f'{method.capitalize()} prediction, {target}; {qualifier}. Stat. errors only; final bin includes overflow.',ha='center',fontsize=11)
            fig.subplots_adjust(top=.90,bottom=.10,left=.075,right=.985)
            save(fig,f'combined_{target}_{method}_summary')
        for cat in ('resolved','merged'):
            for tag in TAGS:
                for ch in CHANNELS:
                    for method,target in [('direct','SR'),('alpha','SR'),('alpha','VR')]:
                        fig=massplot(data,summary,cat,tag,ch,target,method,'inclVR')
                        save(fig,f'{cat}_{tag}_{ch}_{target}_{method}')
                fig=massplot(data,summary,cat,tag,'combined','VR','alpha','heldoutVR')
                save(fig,f'{cat}_{tag}_combined_VR_alpha_heldoutVR')
        for cat in ('resolved','merged'):
            for scope in ('shared',)+TAGS:
                fig,ax=plt.subplots(figsize=(10,7.8));fig.subplots_adjust(left=.15,right=.96,top=.91,bottom=.16)
                ranges=[]
                for cal,color,label in [('inclVR',COLORS['shared'],'CR + VR calibration'),
                                        ('heldoutVR',COLORS['category'],'CR-only calibration')]:
                    root=f'{cat}__SF__{scope}__{cal}__';edges=data[root+'edges'];x=display_edges(edges)
                    values,var=data[root+'values'];error=np.sqrt(var)
                    ax.step(x,step(values),where='post',color=color,lw=1.8,label=label)
                    ax.fill_between(x,step(values-error),step(values+error),step='post',color=color,alpha=.16,linewidth=0)
                    ranges.append(x[-1])
                hep.cms.label('Work in progress',data=True,lumi=59.8,com=13,fontsize=19,ax=ax)
                label='shared b-tag / untagged / VBF-tag' if scope=='shared' else LABELS[scope]
                ax.text(.035,.97,f'2018, {cat}\n{label}\n'+r'$ee+\mu\mu$',transform=ax.transAxes,va='top',fontsize=14)
                ax.legend(loc='upper right',bbox_to_anchor=(.99,.80),fontsize=12)
                ax.set_xlabel(r'$p_{T}^{qq}$ [GeV]' if cat=='resolved' else r'$p_{T}^{J}$ [GeV]',fontsize=20)
                ax.set_ylabel('(Data - TTplusWW - VZ) / DY',fontsize=17)
                ax.axhline(1,color='0.3',lw=.8);ax.set_xlim(0,max(ranges))
                y=ax.get_ylim();ax.set_ylim(min(0,y[0]),y[1]*1.5)
                ax.tick_params(which='both',direction='in',top=True,right=True,labelsize=13)
                fig.text(.5,.035,'Statistical errors only; last pT bin includes the full overflow tail.',ha='center',fontsize=11)
                save(fig,f'{cat}_{scope}_pt_scale_factors')
    shutil.copy2(Path(__file__),output/Path(__file__).name)
    metadata=json.loads((output/'metadata.json').read_text());metadata['figures']=count
    metadata['formats']=['pdf','png'];metadata['style']='mplhep.style.CMS'
    (output/'metadata.json').write_text(json.dumps(metadata,indent=2))
    print(f'Saved {count} individual PNG/PDF pairs and {count}-page PDF bundle',flush=True)


if __name__=='__main__':main()
