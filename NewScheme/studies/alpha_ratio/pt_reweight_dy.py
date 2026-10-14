"""Evaluate category-specific and topology-shared DY pT corrections from saved cells."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.stats import chi2 as chi2_distribution
from pt_reweight_core import (Space,fit_pt,mass_groups,group_edges,predictions,
                             compare,TAGS,CHANNELS,REGIONS)


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--min-effective',type=float,default=25.)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Use a fresh output directory')
    args.output.mkdir(parents=True)
    archive={};rows=[];sfrows=[];binrows=[];binning=[];fitsummary=[]
    with np.load(args.source/'source_histograms.npz') as cached:
        arrays={key:cached[key] for key in cached.files}
    for cat in ('resolved','merged'):
        space=Space(arrays,cat)
        fits={}
        for calibration,regions in [('inclVR',[0,2]),('heldoutVR',[0])]:
            for scope,tags in [('shared',[0,1,2])]+[(t,[i]) for i,t in enumerate(TAGS)]:
                name=scope+'__'+calibration
                fit=fit_pt(space,tags,regions,min_effective=args.min_effective)
                fits[name]=fit
                archive[f'{cat}__SF__{name}__edges']=fit['edges']
                archive[f'{cat}__SF__{name}__values']=np.stack([
                    fit['factor'].value,np.diag(fit['factor'].covariance(space.variance))])
                archive[f'{cat}__SF__{name}__cov']=fit['factor'].covariance(space.variance)
                for i in range(len(fit['groups'])):
                    sfrows.append(dict(category=cat,scope=scope,calibration=calibration,
                        pt_low=float(fit['edges'][i]),pt_high=float(fit['edges'][i+1]),
                        sf=float(fit['factor'].value[i]),
                        sf_stat=float(np.sqrt(fit['factor'].covariance(space.variance)[i,i])),
                        data_minus_background=float(fit['observation'].value[i]),
                        DY=float(fit['dy'].value[i]),effective_Z=float(fit['effective_z'][i]),
                        effective_DY=float(fit['effective_dy'][i])))
                fitsummary.append(dict(category=cat,scope=scope,calibration=calibration,
                    pt_bins=len(fit['groups']),pt_edges=[None if np.isinf(x) else float(x) for x in fit['edges']]))
        for tagindex,tag in enumerate(TAGS):
            relevant=[fits[scope+'__'+cal] for scope in (tag,'shared') for cal in ('inclVR','heldoutVR')]
            groups=mass_groups(space,tagindex,relevant)
            edges=group_edges(space.mass_edges,groups)
            archive[f'{cat}__{tag}__mass_edges']=edges
            binning.append(dict(category=cat,tag=tag,edges=[None if np.isinf(x) else float(x) for x in edges]))
            for channel,chs in [('ee',[0]),('mumu',[1]),('combined',[0,1])]:
                for target,region in [('SR',1),('VR',2)]:
                    datacounts=np.array([space.data[np.ix_([tagindex],chs,[region],mg,np.arange(space.shape[-1]))].sum() for mg in groups])
                    nominaldy=space.linear(1,[tagindex],chs,[region],groups)
                    var=np.diag(nominaldy.covariance(space.variance))
                    effective=np.divide(nominaldy.value**2,var,out=np.zeros_like(var),where=var>0)
                    use=(datacounts>=5)&(effective>=5)
                    for scope in ('category','shared'):
                        for calibration in ('inclVR','heldoutVR'):
                            fit=fits[(tag if scope=='category' else 'shared')+'__'+calibration]
                            result=predictions(space,tagindex,chs,region,groups,fit)
                            if target=='SR':
                                inclusive=predictions(space,tagindex,chs,region,groups,fit,cr_regions=(0,2))
                                result['alpha_inclusiveCR_nominal']=inclusive['alpha_nominal']
                                result['alpha_inclusiveCR_reweighted']=inclusive['alpha_reweighted']
                            prefix='__'.join((cat,tag,channel,target,scope,calibration))+'__'
                            archive[prefix+'edges']=edges
                            archive[prefix+'chi2_mask']=use
                            for name,prediction in result.items():
                                cov=prediction.covariance(space.variance)
                                archive[prefix+name]=np.stack([prediction.value,np.diag(cov)])
                                archive[prefix+name+'_cov']=cov
                                if name.startswith('transfer') or name=='observed':continue
                                residual=result['observed']-prediction
                                archive[prefix+name+'_residual_cov']=residual.covariance(space.variance)
                                ratio=np.full((2,len(groups)),np.nan)
                                good=prediction.value>0
                                if good.any():
                                    from pt_reweight_core import Estimate
                                    numerator=Estimate(result['observed'].value[good],result['observed'].jac[good])
                                    denominator=Estimate(prediction.value[good],prediction.jac[good])
                                    rr=numerator/denominator
                                    ratio[0,good]=rr.value
                                    ratio[1,good]=np.diag(rr.covariance(space.variance))
                                archive[prefix+name+'_closure_bins']=ratio
                                row=dict(category=cat,tag=tag,channel=channel,target=target,
                                         scope=scope,calibration=calibration,prediction=name)
                                row.update(compare(result['observed'],prediction,space.variance,valid=use))
                                row['chi2_pvalue']=float(chi2_distribution.sf(row['chi2'],row['rank'])) if row['rank'] else np.nan
                                row['calibration_includes_target']=(target=='VR' and calibration=='inclVR')
                                rows.append(row)
                                for i in range(len(groups)):
                                    binrows.append(dict(category=cat,tag=tag,channel=channel,target=target,
                                        scope=scope,calibration=calibration,prediction=name,
                                        low=float(edges[i]),high=float(edges[i+1]),
                                        observation=float(result['observed'].value[i]),
                                        prediction_value=float(prediction.value[i]),variance=float(cov[i,i]),
                                        residual_variance=float(archive[prefix+name+'_residual_cov'][i,i]),
                                        closure=float(ratio[0,i]),closure_variance=float(ratio[1,i]),
                                        chi2_used=bool(use[i])))
                    print(f'{cat} {tag} {channel} {target}: {len(groups)} bins',flush=True)
        del space
    write_csv(args.output/'summary.csv',rows)
    write_csv(args.output/'pt_factors.csv',sfrows)
    write_csv(args.output/'bins.csv',binrows)
    np.savez_compressed(args.output/'results.npz',**archive)
    metadata=dict(date='2026-10-07',source=str(args.source.resolve()),
        source_sha256=hashlib.sha256((args.source/'source_histograms.npz').read_bytes()).hexdigest(),
        ZH_subtracted=False,min_effective_pt=args.min_effective,
        calibrations={'inclVR':'CR+VR (all non-SR hadronic masses 40-180)',
                      'heldoutVR':'CR only (both SR and VR held out)'},
        alpha_source='Primary SR uses CR+VR (all non-SR); primary VR uses CR only (disjoint target). Also save SR CR-only sensitivity.',
        correction='Shared ee+mumu fit in each category or topology, absolute binned (Data-TTplusWW-VZ)/DY',
        uncertainties='Joint first-order Jacobian covariance from independent signed-Z and DY cells; includes SF calibration, MC, CR, shared channels and VR reuse',
        missing_systematics='Lumi, tags, subtraction normalization/shape, DY modeling and normalization-induced MC correlations; adaptive bin choice treated as fixed',
        chi2='Gaussian full-covariance diagnostic, excludes bins with target data count<5 or nominal DY Neff<5; not a final likelihood test',
        mass_binning='Merge adjacent bins until nominal and reweighted DY positive separately in channels/regions; no SR/VR target data used for mass-bin choice',
        fit_binning='Merge adjacent pT bins based only on calibration signed yields and effective counts, no target-SR observations',
        fits=fitsummary,mass_binning_edges=binning,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output/'metadata.json').write_text(json.dumps(metadata,indent=2))
    shutil.copy2(args.source/'source_audit.json',args.output/'source_audit.json')
    for filename in ('pt_reweight_scan.py','pt_reweight_core.py','pt_reweight_dy.py','test_pt_reweight.py'):
        shutil.copy2(Path(__file__).parent/filename,args.output/filename)
    print('Analysis complete: '+str(args.output),flush=True)


if __name__=='__main__':main()
