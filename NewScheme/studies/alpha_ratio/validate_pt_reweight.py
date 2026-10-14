"""Validate the source scan against existing reductions and test shared-tag calibration."""
import argparse
import csv
import json
import os
from pathlib import Path
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import numpy as np
from scipy.stats import chi2 as chi2_distribution
from pt_reweight_core import Space,fit_pt,compare,TAGS,Estimate


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--study',type=Path,required=True)
    p.add_argument('--baseline',type=Path,required=True)
    args=p.parse_args()
    with np.load(args.source/'source_histograms.npz') as f:arrays={key:f[key] for key in f.files}
    checked=0;relative=0.
    with np.load(args.baseline/'histograms.npz') as old:
        for cat in ('resolved','merged'):
            for tag in TAGS:
                for channel,ch in enumerate(('ee','mumu')):
                    for ri,region in enumerate(('CR','SR','VR')):
                        prefix=f'{cat}__{tag}__{ch}__'+('SR' if region=='CR' else region)+'__'
                        for proc in ('Data','DY','TTplusWW','VZ'):
                            key=prefix+('CR_' if region=='CR' else 'target_')+proc
                            expected=old[key].sum(axis=1)
                            actual=arrays[f'{cat}__{tag}__{proc}'][channel,ri].sum(axis=(0,1))
                            np.testing.assert_allclose(actual,expected,rtol=2e-6,atol=2e-5,
                                err_msg=f'{cat} {tag} {ch} {region} {proc}')
                            relative=max(relative,float(np.max(np.abs(actual-expected)/np.maximum(np.abs(expected),1))))
                            checked+=1
    rows=[];global_rows=[]
    for cat in ('resolved','merged'):
        space=Space(arrays,cat)
        for calibration,regions in [('inclVR',[0,2]),('heldoutVR',[0])]:
            fit=fit_pt(space,[0,1,2],regions)
            total_w=space.weighted([0,1,2],[0,1],regions,[np.arange(space.shape[-2])],fit['groups'],fit['factor'])
            total_z=space.linear(0,[0,1,2],[0,1],regions)
            np.testing.assert_allclose(total_w.value,total_z.value,rtol=1e-10,atol=1e-8)
            np.testing.assert_allclose(total_w.jac,total_z.jac,rtol=1e-9,atol=1e-8)
            residuals=[];valids=[]
            for ti,tag in enumerate(TAGS):
                z=space.linear(0,[ti],[0,1],regions,ptgroups=fit['groups'])
                dy=space.linear(1,[ti],[0,1],regions,ptgroups=fit['groups'])
                predicted=dy*fit['factor']
                vz=np.diag(z.covariance(space.variance));vd=np.diag(dy.covariance(space.variance))
                effective_z=np.divide(z.value**2,vz,out=np.zeros_like(vz),where=vz>0)
                effective_dy=np.divide(dy.value**2,vd,out=np.zeros_like(vd),where=vd>0)
                use=(effective_z>=5)&(effective_dy>=5)&(dy.value>0)
                row=dict(category=cat,tag=tag,calibration=calibration)
                row.update(compare(z,predicted,space.variance,use))
                row['chi2_pvalue']=float(chi2_distribution.sf(row['chi2'],row['rank'])) if row['rank'] else np.nan
                rows.append(row);residuals.append(z-predicted);valids.append(use)
            residual=Estimate(np.concatenate([r.value for r in residuals]),np.concatenate([r.jac for r in residuals]))
            zero=Estimate(np.zeros_like(residual.value),np.zeros_like(residual.jac))
            check=compare(residual,zero,space.variance,np.concatenate(valids))
            global_rows.append(dict(category=cat,calibration=calibration,chi2=check['chi2'],rank=check['rank'],
                pvalue=float(chi2_distribution.sf(check['chi2'],check['rank'])) if check['rank'] else None))
    with (args.study/'shared_tag_compatibility.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    report=dict(source_vs_reduced_comparisons=checked,maximum_relative_difference=relative,
                source_vs_reduced='passed',calibration_closure_mean_and_jacobian='passed',
                common_tag_correction_tests=global_rows)
    (args.study/'validation.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
