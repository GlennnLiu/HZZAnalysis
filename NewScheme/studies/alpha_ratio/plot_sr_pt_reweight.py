"""Supplement the direct-SR study with the hadronic-pT distribution itself."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import mplhep as hep
from scipy.stats import chi2 as chi2_distribution
from pt_reweight_core import Space,Estimate,fit_pt,compare,TAGS
from plot_pt_reweight import massplot,prefix,display_edges


def aggregate(estimate,groups):
    return Estimate(np.array([estimate.value[g].sum() for g in groups]),
                    np.stack([estimate.jac[g].sum(axis=0) for g in groups]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--study',type=Path,required=True)
    args=parser.parse_args()
    if (args.study/'all_SR_pt_reweight.pdf').exists():raise FileExistsError('Supplement already exists')
    with np.load(args.source/'source_histograms.npz') as f:arrays={k:f[k] for k in f.files}
    with (args.study/'summary.csv').open() as f:mainrows=list(csv.DictReader(f))
    summary={tuple(r[k] for k in ('category','tag','channel','target','scope','calibration','prediction')):r for r in mainrows}
    plt.style.use(hep.style.CMS)
    plt.rcParams.update({'font.size':14,'axes.labelpad':7,'pdf.fonttype':42,'axes.linewidth':1.1})
    saved={};rows=[]
    with PdfPages(args.study/'all_SR_pt_reweight.pdf') as bundle:
        for cat in ('resolved','merged'):
            space=Space(arrays,cat)
            fits={scope:fit_pt(space,([0,1,2] if scope=='shared' else [TAGS.index(scope)]),[0,2])
                  for scope in ('shared',)+TAGS}
            pg=fits['shared']['groups'];edges=fits['shared']['edges']
            finegroups=[np.array([i]) for i in range(space.shape[-1])]
            for ti,tag in enumerate(TAGS):
                for ch,channels in [('ee',[0]),('mumu',[1]),('combined',[0,1])]:
                    obs=space.linear(0,[ti],channels,[1],ptgroups=pg)
                    fine_dy=space.linear(1,[ti],channels,[1],ptgroups=finegroups)
                    nominal=aggregate(fine_dy,pg)
                    variance=np.diag(nominal.covariance(space.variance))
                    eff=np.divide(nominal.value**2,variance,out=np.zeros_like(variance),where=variance>0)
                    counts=np.array([space.data[np.ix_([ti],channels,[1],np.arange(space.shape[-2]),g)].sum() for g in pg])
                    use=(counts>=5)&(eff>=5)
                    for scope in ('category','shared'):
                        fit=fits[tag if scope=='category' else 'shared']
                        lookup=np.zeros(space.shape[-1],dtype=int)
                        for i,g in enumerate(fit['groups']):lookup[g]=i
                        expanded=Estimate(fit['factor'].value[lookup],fit['factor'].jac[lookup])
                        weighted=aggregate(fine_dy*expanded,pg)
                        norm=fit['observation'].total()/fit['dy'].total()
                        norm=Estimate(np.repeat(norm.value,len(pg)),np.repeat(norm.jac,len(pg),axis=0))
                        items={'observed':obs,'dy_nominal':nominal,'dy_reweighted':weighted,'dy_normonly':nominal*norm}
                        root=prefix(cat,tag,ch,'SR',scope,'inclVR')
                        saved[root+'edges']=edges;saved[root+'chi2_mask']=use
                        for key,item in items.items():
                            saved[root+key]=np.stack([item.value,np.diag(item.covariance(space.variance))])
                            saved[root+key+'_cov']=item.covariance(space.variance)
                            if key=='observed':continue
                            result=compare(obs,item,space.variance,use)
                            expected=summary[(cat,tag,ch,'SR',scope,'inclVR',key)]
                            np.testing.assert_allclose(result['predicted'],float(expected['predicted']),rtol=1e-10,atol=1e-7)
                            np.testing.assert_allclose(result['predicted_stat'],float(expected['predicted_stat']),rtol=1e-9,atol=1e-7)
                            result['chi2_pvalue']=float(chi2_distribution.sf(result['chi2'],result['rank'])) if result['rank'] else np.nan
                            rows.append(dict(category=cat,tag=tag,channel=ch,scope=scope,prediction=key,**result))
                            good=item.value>0
                            ratio=np.full((2,len(pg)),np.nan)
                            if good.any():
                                rr=Estimate(obs.value[good],obs.jac[good])/Estimate(item.value[good],item.jac[good])
                                ratio[0,good]=rr.value;ratio[1,good]=np.diag(rr.covariance(space.variance))
                            saved[root+key+'_closure_bins']=ratio
                    fig=massplot(saved,summary,cat,tag,ch,'SR','direct','inclVR')
                    lower=fig.axes[1];x=display_edges(edges);center=(x[:-1]+x[1:])/2
                    ticks=[v for v in (0,200,300,400,500,650,850) if v<x[-2]]+[center[-1]]
                    lower.set_xticks(ticks);lower.set_xticklabels([f'{v:g}' for v in ticks[:-1]]+[rf'$\geq {x[-2]:g}$'])
                    lower.set_xlabel(r'$p_{T}^{qq}$ [GeV]' if cat=='resolved' else r'$p_{T}^{J}$ [GeV]',fontsize=20)
                    stem=f'{cat}_{tag}_{ch}_SR_direct_pt'
                    fig.savefig(args.study/(stem+'.pdf'),bbox_inches='tight')
                    fig.savefig(args.study/(stem+'.png'),dpi=180,bbox_inches='tight')
                    bundle.savefig(fig,bbox_inches='tight');plt.close(fig)
    with (args.study/'sr_pt_summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    np.savez_compressed(args.study/'sr_pt_distributions.npz',**saved)
    shutil.copy2(Path(__file__),args.study/Path(__file__).name)
    metadata=json.loads((args.study/'metadata.json').read_text())
    metadata['supplementary_SR_pt_figures']=18;metadata['total_individual_figures']=89
    (args.study/'metadata.json').write_text(json.dumps(metadata,indent=2))
    print('Saved 18 SR pT PNG/PDF pairs; all weighted yields/errors agree with search-mass projections.',flush=True)


if __name__=='__main__':main()
