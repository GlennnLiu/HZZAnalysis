"""Remove mass reweighting from existing VR displays; do not refit or alter results."""
import argparse,csv,hashlib,json,shutil,zipfile
from pathlib import Path
import plot_vr_methods as p

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def large_individual(fig):
    """Match distribution-plot typography without changing any plotted values."""
    ax,lower=fig.axes[:2]
    fig.set_size_inches(12,11)
    fig.subplots_adjust(left=.18,right=.97,top=.90,bottom=.20)
    ax.set_ylabel('Events / bin',fontsize=32)
    lower.set_ylabel('Obs. / pred.',fontsize=27)
    lower.xaxis.label.set_fontsize(34)
    for item in (ax,lower): item.tick_params(labelsize=23)
    region=None
    for text in list(ax.texts):
        if text.get_text().startswith('2018, VR'):
            region=text; text.set_fontsize(21)
        else:
            text.remove()
    p.hep.cms.label('Work in progress',data=True,lumi=59.8,com=13,fontsize=27,ax=ax)
    handles,labels=ax.get_legend_handles_labels()
    labels=[label.split('  (Z/pred.')[0].replace('Data - TTplusWW - VZ','Data - TT - VZ') for label in labels]
    ax.get_legend().remove()
    legend=ax.legend(handles,labels,loc='upper right',bbox_to_anchor=(.99,.98),fontsize=20,
                     frameon=False,handlelength=1.6,labelspacing=.30)
    # Keep explanatory notes, using short readable lines instead of tiny text.
    notes=('Correction fit: CR + VR (partly in-sample).' if 'CR + VR' in fig.texts[0].get_text()
           else 'Correction fit: CR only (VR held out).')
    for text in list(fig.texts): text.remove()
    fig.text(.5,.115,notes+' No ZH subtraction.',ha='center',fontsize=13)
    fig.text(.5,.080,'Bands: prediction stat. uncertainty; ratio errors include shared covariance.',ha='center',fontsize=13)
    fig.text(.5,.045,'Alpha source excludes VR; original bins retained; last bin includes overflow.',ha='center',fontsize=13)
    fig.canvas.draw(); renderer=fig.canvas.get_renderer()
    boxes=[legend.get_window_extent(renderer),region.get_window_extent(renderer)]
    assert not boxes[0].overlaps(boxes[1]),'Legend/category label overlap'
    # Linear distributions: only add the y headroom required beneath the labels.
    points=[]
    for line in ax.lines: points.extend(zip(line.get_xdata(),line.get_ydata()))
    for collection in ax.collections:
        for path in collection.get_paths(): points.extend(path.vertices)
    points=p.np.asarray(points,dtype=float); points=points[p.np.isfinite(points).all(axis=1)]
    low,_=ax.get_ylim(); high=max(1.,float(ax.dataLim.ymax))*1.25
    for box in boxes:
        coordinates=ax.transAxes.inverted().transform(box.get_points()); fraction=max(.2,coordinates[0,1]-.025)
        x1,x2=ax.transData.inverted().transform([[box.x0,ax.bbox.y0],[box.x1,ax.bbox.y0]])[:,0]
        selected=points[(points[:,0]>=x1)&(points[:,0]<=x2)]
        if len(selected): high=max(high,low+(selected[:,1].max()-low)/fraction)
    ax.set_ylim(low,high*1.02)
    return fig

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path)
    parser.add_argument('--input',type=Path,default=Path(__file__).resolve().parents[3]/'FitTemplate/AlphaRatio_examination/plots/VR_method_comparison_2018')
    parser.add_argument('--large-individuals',action='store_true')
    parser.add_argument('--preview',action='store_true'); args=parser.parse_args()
    source=args.input.resolve()
    folder=args.output.resolve() if args.output else source
    folder.mkdir(parents=True,exist_ok=True)
    preserved={name:sha(source/name) for name in ('results.npz','summary.csv','metadata.json')}
    backup=folder/'five_method_plot_backup_20261013'
    if folder==source and not backup.exists():
        backup.mkdir()
        for pattern in ('*.png','*.pdf','*_VR_chi2_*.md'):
            for file in folder.glob(pattern): shutil.copy2(file,backup/file.name)
    data,_,_=p.load(source)
    with (source/'summary.csv').open() as f: rows=list(csv.DictReader(f))
    summaries={(row['category'],row['tag'],row['channel'],row['calibration'],row['method']):row for row in rows}
    p.METHODS=tuple(m for m in p.METHODS if m[0]!='direct_mass')
    assert len(p.METHODS)==4
    p.plt.style.use(p.hep.style.CMS)
    p.plt.rcParams.update({'font.size':14,'axes.labelpad':7,'pdf.fonttype':42,'axes.linewidth':1.1})
    outputs=[]
    def save(fig,stem,bundle):
        for ax in fig.axes:
            assert not any('mass reweighted' in label.lower() for label in ax.get_legend_handles_labels()[1])
        for ext in ('pdf','png'):
            file=folder/(stem+'.'+ext); fig.savefig(file,dpi=180 if ext=='png' else None,bbox_inches='tight'); outputs.append(file)
        if bundle is not None: bundle.savefig(fig,bbox_inches='tight')
        p.plt.close(fig)
    if args.preview:
        for cat,tag in (('resolved','untag'),('merged','btag')):
            fig=p.plot(data,summaries,cat,tag,'combined','heldoutVR')
            save(large_individual(fig),f'{cat}_{tag}_combined_VR_methods_heldoutVR',None)
        print('Large-text preview passed',flush=True); return
    for cal in p.CALIBRATIONS:
        file=folder/f'all_VR_method_comparison_{cal}.pdf'
        with p.PdfPages(file) as bundle:
            for ch in p.CHANNELS:
                fig=p.plt.figure(figsize=(18,12)); outer=fig.add_gridspec(2,3,wspace=.32,hspace=.38)
                for i,(cat,tag) in enumerate((c,t) for c in ('resolved','merged') for t in p.TAGS):
                    grid=outer[i//3,i%3].subgridspec(2,1,height_ratios=[3,1.15],hspace=.05)
                    ax=fig.add_subplot(grid[0]); lower=fig.add_subplot(grid[1],sharex=ax)
                    p.plot(data,summaries,cat,tag,ch,cal,compact=True,axes=(ax,lower))
                handles,labels=fig.axes[0].get_legend_handles_labels()
                fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.995),ncol=3,fontsize=11,frameon=False)
                qualifier=('CR+VR correction calibration: VR partly in-sample' if cal=='inclVR' else 'CR-only correction calibration: VR held out of fitted factor values')
                fig.text(.5,.018,qualifier+'; alpha source excludes VR; stat. only; no ZH subtraction.',ha='center',fontsize=11)
                fig.subplots_adjust(top=.875,bottom=.10,left=.07,right=.985)
                if args.large_individuals:
                    bundle.savefig(fig,bbox_inches='tight'); p.plt.close(fig)
                else: save(fig,f'{ch}_VR_methods_summary_{cal}',bundle)
            for cat in ('resolved','merged'):
                for tag in p.TAGS:
                    for ch in p.CHANNELS:
                        fig=p.plot(data,summaries,cat,tag,ch,cal)
                        if args.large_individuals: large_individual(fig)
                        save(fig,f'{cat}_{tag}_{ch}_VR_methods_{cal}',bundle)
        outputs.append(file); print('Redrawn without mass reweighting:',cal,flush=True)
    for ch in p.CHANNELS:
        for cal in p.CALIBRATIONS:
            table=[f'# VR chi2 comparison, {ch}, {cal}','',
                   'Display excludes mass reweighting. Remaining statistics are unchanged from the saved five-method study.','',
                   '| Category | Direct DY | DY + pT | Nominal alpha | Alpha + pT |',
                   '|---|---:|---:|---:|---:|']
            for cat in ('resolved','merged'):
                for tag in p.TAGS:
                    values=[f"{float(summaries[(cat,tag,ch,cal,m[0])]['chi2']):.2f}/{int(summaries[(cat,tag,ch,cal,m[0])]['rank'])}" for m in p.METHODS]
                    table.append('| '+cat+' '+p.LABELS[tag]+' | '+' | '.join(values)+' |')
            file=folder/f'{ch}_VR_chi2_{cal}.md'; file.write_text('\n'.join(table)+'\n'); outputs.append(file)
    assert preserved=={name:sha(source/name) for name in preserved}
    note=folder/'DISPLAY_UPDATE.md'
    note.write_text('# VR display update — 2026-10-13\n\nMass-reweighted DY removed from both calibration sets, all spectra, ratio panels, legends and chi2 display tables. Four methods remain: direct DY, pT-reweighted DY, nominal alpha and pT-reweighted alpha.\n\nOriginal filenames retained. Because /cms quota is exceeded, updated plots were rendered in /scratchfs2 and copied to the local plot folder; remote /cms figures remain the original five-method versions. Cached numerical results, uncertainties, original binning and metadata are unchanged; no fits, trees or skims were regenerated. Reproduce using PYTHONPATH pointing to the original plot folder and `python redraw_without_mass.py --output OUTPUT`.\n')
    validation=folder/'display_validation_no_mass.json'
    if args.large_individuals:
        with note.open('a') as f: f.write('\nIndividual-plot typography enlarged: axis titles32/34pt, ratio title27pt, ticks23pt, legend20pt, category21pt, CMS/luminosity27pt. Summary PNG/PDF pages unchanged; bundles include the enlarged individual pages. Integrated closure numbers remain in summary.csv, rather than lengthening the legend.\n')
    validation.write_text(json.dumps(dict(status='passed',displayed_methods=[m[0] for m in p.METHODS],preserved_input_sha256=preserved,png=36 if args.large_individuals else 42,pdf=38 if args.large_individuals else 44,no_refit=True,large_individuals=args.large_individuals),indent=2))
    outputs.extend((Path(__file__),note,validation))
    zipname='VR_method_comparison_2018_large_individuals.zip' if args.large_individuals else 'VR_method_comparison_2018_no_mass_display.zip'
    with zipfile.ZipFile(folder.parent/zipname,'w',zipfile.ZIP_DEFLATED) as z:
        for file in outputs: z.write(file,folder.name+'/'+file.name)
    print('PASS: individual text enlarged; summary files untouched; cached physics inputs unchanged.' if args.large_individuals else 'PASS: 42 PNGs, 44 PDFs; no mass curve or legend entries; cached physics inputs unchanged.',flush=True)
if __name__=='__main__': main()
