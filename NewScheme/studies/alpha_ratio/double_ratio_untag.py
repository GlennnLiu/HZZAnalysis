"""Double ratios and sideband-only DY transfer corrections; 2018 untagged."""
import argparse
import csv
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(name, '1')
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

QEDGES = np.array([40,55,70,85,105,120,135,150,165,180.])
TRAIN = np.array([0,1,4,5,8])
TARGETS = {'SR': np.array([2,3]), 'VR': np.array([6,7])}
MODELS = ['separable', 'cubic_separable', 'linear_interaction', 'quadratic_interaction', 'local_interpolation']


def read_cells(base, broad, cat, channel, xe):
    components = {}
    for proc in ['Data', 'TTplusWW', 'VZ', 'DY']:
        parts = []
        for reg, directory in [('CR',broad),('SR',base)]:
            path=directory/f'{proc}_{channel}_{cat}_untag_2018_{reg}.root'
            with uproot.open(path) as f:
                parts.append(f['tree'].arrays(['mass','mass2q','weight'],library='np'))
        m=np.concatenate([p['mass'] for p in parts])
        q=np.concatenate([p['mass2q'] for p in parts])
        w=np.concatenate([p['weight'] for p in parts])
        components[proc]=(np.histogram2d(m,q,bins=[xe,QEDGES],weights=w)[0],
                          np.histogram2d(m,q,bins=[xe,QEDGES],weights=w*w)[0])
    z=components['Data'][0]-components['TTplusWW'][0]-components['VZ'][0]
    vz=components['Data'][1]+components['TTplusWW'][1]+components['VZ'][1]
    d,vd=components['DY']
    return z,vz,d,vd,components


def design(x,q,model):
    u=(q-100)/60
    v=np.log(x/500)
    fields=[np.ones_like(u),u,u*u,v,v*v]
    if model == 'cubic_separable': fields.append(u*u*u)
    if model in ['linear_interaction','quadratic_interaction']: fields.append(u*v)
    if model == 'quadratic_interaction': fields.append(u*u*v)
    return np.stack(fields,axis=-1)


def fit_surface(z,vz,d,vd,xe,model):
    xc=(xe[:-1]+xe[1:])/2
    qc=(QEDGES[:-1]+QEDGES[1:])/2
    if model == 'local_interpolation':
        # Fit each sideband slice versus total mass, then interpolate log k
        # locally in hadronic mass. SR and VR cells never enter these fits.
        v=np.log(xc/500)
        anchors=[]; totalx2=0.; totalndof=0
        for j in TRAIN:
            good=(z[:,j]>0)&(d[:,j]>0)&(z[:,j]**2/np.maximum(vz[:,j],1e-20)>=4)&(d[:,j]**2/np.maximum(vd[:,j],1e-20)>=4)
            degree=min(2,int(good.sum())-2)
            if degree<0: raise ValueError('Too few local anchor cells')
            basis=np.stack([v**i for i in range(degree+1)],axis=-1)
            variance=vz[good,j]/z[good,j]**2+vd[good,j]/d[good,j]**2
            a=basis[good]/np.sqrt(variance)[:,None]
            y=np.log(z[good,j]/d[good,j])/np.sqrt(variance)
            coeff=np.linalg.lstsq(a,y,rcond=None)[0]
            anchors.append(basis@coeff)
            residual=y-a@coeff
            totalx2+=float(residual@residual);totalndof+=int(good.sum()-len(coeff))
        logk=np.stack([np.interp(qc,qc[TRAIN],np.array(anchors)[:,i]) for i in range(len(xc))])
        matrix=np.eye(logk.size).reshape(*logk.shape,logk.size)
        return logk.ravel(),np.eye(logk.size),totalx2,totalndof,matrix
    xx,qq=np.meshgrid(xc,qc,indexing='ij')
    matrix=design(xx,qq,model)
    train=np.zeros_like(z,dtype=bool); train[:,TRAIN]=True
    valid=train&(z>0)&(d>0)&(vz>0)
    # Avoid taking logs of cells dominated by signed-weight cancellation or
    # very low effective statistics. No clipping negative yields into positives.
    valid &= (z*z/vz>=4)&(d*d/np.maximum(vd,1e-20)>=4)
    if valid.sum()<=matrix.shape[-1]: raise ValueError('Insufficient sideband cells')
    variance=vz[valid]/z[valid]**2+vd[valid]/d[valid]**2
    a=matrix[valid]/np.sqrt(variance)[:,None]
    y=np.log(z[valid]/d[valid])/np.sqrt(variance)
    coeff=np.linalg.lstsq(a,y,rcond=None)[0]
    cov=np.linalg.inv(a.T@a)
    residual=y-a@coeff
    return coeff,cov,float(residual@residual),int(valid.sum()-len(coeff)),matrix


def predictions(z,d,k):
    c=z[:,TRAIN].sum(axis=1)
    dc=d[:,TRAIN].sum(axis=1)
    kc=(d[:,TRAIN]*k[:,TRAIN]).sum(axis=1)/dc
    out={}
    for target,indices in TARGETS.items():
        dt=d[:,indices].sum(axis=1)
        kt=(d[:,indices]*k[:,indices]).sum(axis=1)/dt
        out[target]=c*dt/dc*kt/kc
    return out


def savefig(fig,output,stem):
    for ext in ['png','pdf']: fig.savefig(output/f'{stem}.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base',type=Path,required=True)
    p.add_argument('--cr-base',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--bootstrap',type=int,default=500)
    args=p.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(20261005)
    cells,drrows,summary,fitrows,drtests,predrows,regionrows=[],[],[],[],[],[],[]
    study={}
    for cat in ['resolved','merged']:
        xe=np.array([300,400,500,650,850,1100,1600,3500.]) if cat=='resolved' else np.array([300,500,850,3500.])
        channels={}
        for channel in ['ee','mumu']:
            channels[channel]=read_cells(args.base,args.cr_base,cat,channel,xe)
        for channel in ['ee','mumu','combined']:
            if channel=='combined':
                z,vz,d,vd=[channels['ee'][i]+channels['mumu'][i] for i in range(4)]
            else: z,vz,d,vd=channels[channel][:4]
            nz=len(xe)-1
            xc=(xe[:-1]+xe[1:])/2
            qc=(QEDGES[:-1]+QEDGES[1:])/2
            r=z/d
            lrvar=vz/np.maximum(z*z,1e-20)+vd/np.maximum(d*d,1e-20)
            # Four independent cells: x_i/q_j, x_0/q_j, x_i/q_low, x_0/q_low.
            # For x_0 or q_low, numerator and denominator coincide exactly.
            dr=r/r[0:1,:]/(r[:,0:1]/r[0,0])
            vdr=dr*dr*(lrvar+lrvar[0:1,:]+lrvar[:,0:1]+lrvar[0,0])
            dr[0,:]=1; dr[:,0]=1; vdr[0,:]=0; vdr[:,0]=0
            dr[(z<=0)|(d<=0)]=np.nan
            for label,js in [('sidebands',[1,4,5,8]),('SR',[2,3]),('VR',[6,7])]:
                contrasts=[];values=[]
                for i in range(1,nz):
                    for j in js:
                        indices=[(i,j),(0,j),(i,0),(0,0)]
                        if not all(z[a,b]>0 and d[a,b]>0 and z[a,b]**2/max(vz[a,b],1e-20)>=4 and d[a,b]**2/max(vd[a,b],1e-20)>=4 for a,b in indices):continue
                        contrast=np.zeros(z.size)
                        for (a,b),sign in zip(indices,[1,-1,-1,1]):contrast[a*z.shape[1]+b]+=sign
                        contrasts.append(contrast);values.append(np.log(dr[i,j]))
                if contrasts:
                    matrixc=np.array(contrasts)
                    covariance=(matrixc*lrvar.ravel())@matrixc.T
                    value=np.array(values)
                    drtests.append(dict(category=cat,channel=channel,region=label,stat_chi2=float(value@np.linalg.pinv(covariance)@value),rank=int(np.linalg.matrix_rank(covariance)),contrasts=len(value)))
            for i in range(nz):
                for j in range(len(qc)):
                    cells.append(dict(category=cat,channel=channel,xlow=xe[i],xhigh=xe[i+1],qlow=QEDGES[j],qhigh=QEDGES[j+1],signed=z[i,j],signed_variance=vz[i,j],dy=d[i,j],dy_variance=vd[i,j],data_mc=r[i,j]))
                    drrows.append(dict(category=cat,channel=channel,xlow=xe[i],xhigh=xe[i+1],qlow=QEDGES[j],qhigh=QEDGES[j+1],double_ratio=dr[i,j],stat=np.sqrt(vdr[i,j])))
            fig,ax=plt.subplots(figsize=(10,5))
            display=np.ma.masked_invalid(dr.T)
            im=ax.imshow(display,origin='lower',aspect='auto',vmin=.7,vmax=1.3,cmap='coolwarm')
            ax.set_xticks(np.arange(nz)); ax.set_xticklabels([f'{a:g}-{b:g}' for a,b in zip(xe[:-1],xe[1:])],rotation=30)
            ax.set_yticks(np.arange(len(qc))); ax.set_yticklabels([f'{a:g}-{b:g}' for a,b in zip(QEDGES[:-1],QEDGES[1:])])
            for i in range(nz):
                for j in range(len(qc)):
                    if np.isfinite(dr[i,j]): ax.text(i,j,f'{dr[i,j]:.2f}\n+/-{np.sqrt(vdr[i,j]):.2f}',ha='center',va='center',fontsize=7)
            ax.set_xlabel('m(2l2q) interval [GeV]'); ax.set_ylabel('m(qq) interval [GeV]')
            ax.set_title(f'2018 {cat} untag, {channel}: four-region double ratio\nReferences: first total-mass bin and 40-55 GeV hadronic mass')
            fig.colorbar(im,ax=ax,label='Double ratio (shared reference errors correlated)')
            savefig(fig,args.output,f'{cat}_{channel}_double_ratio')
            nominal=predictions(z,d,np.ones_like(z))
            # A direct SR/VR double ratio uses no total-mass control interval.
            # It tests whether the measured VR alpha correction transfers to SR.
            sz=z[:,TARGETS['SR']].sum(axis=1);vs=vz[:,TARGETS['SR']].sum(axis=1)
            vzjet=z[:,TARGETS['VR']].sum(axis=1);vv=vz[:,TARGETS['VR']].sum(axis=1)
            sd=d[:,TARGETS['SR']].sum(axis=1);vsd=vd[:,TARGETS['SR']].sum(axis=1)
            vdy=d[:,TARGETS['VR']].sum(axis=1);vvd=vd[:,TARGETS['VR']].sum(axis=1)
            cross=(sz/sd)/(vzjet/vdy)
            crossvar=cross**2*(vs/np.maximum(sz**2,1e-20)+vsd/np.maximum(sd**2,1e-20)+vv/np.maximum(vzjet**2,1e-20)+vvd/np.maximum(vdy**2,1e-20))
            fig,ax=plt.subplots(figsize=(8,4))
            ax.errorbar(xc,cross,np.sqrt(crossvar),fmt='o',color='black')
            ax.axhline(1,color='gray',ls='--');ax.set_ylim(.5,1.8)
            ax.set_xlabel('m(2l2q) [GeV]');ax.set_ylabel('(signed SR / signed VR) / (DY SR / DY VR)')
            ax.set_title(f'2018 {cat} untag, {channel}: SR/VR double ratio\nDiagnostic only: SR observations do not calibrate the correction')
            savefig(fig,args.output,f'{cat}_{channel}_SR_VR_double_ratio')
            for i in range(nz):regionrows.append(dict(category=cat,channel=channel,xlow=xe[i],xhigh=xe[i+1],double_ratio=cross[i],stat=np.sqrt(crossvar[i])))
            vc=vz[:,TRAIN].sum(axis=1); c=z[:,TRAIN].sum(axis=1)
            dc=d[:,TRAIN].sum(axis=1); vdc=vd[:,TRAIN].sum(axis=1)
            coeffs={}; boot={model:{target:[] for target in TARGETS} for model in MODELS}
            for model in MODELS:
                try: coef,cov,fitx2,ndof,matrix=fit_surface(z,vz,d,vd,xe,model)
                except (ValueError,np.linalg.LinAlgError): continue
                coeffs[model]=(coef,matrix)
                fitrows.append(dict(category=cat,channel=channel,model=model,train_chi2=fitx2,train_ndof=ndof))
            for b in range(args.bootstrap):
                zb=rng.normal(z,np.sqrt(vz)); db=rng.normal(d,np.sqrt(vd))
                if np.any(db[:,TRAIN].sum(axis=1)<=0): continue
                if any(np.any(db[:,ix].sum(axis=1)<=0) for ix in TARGETS.values()): continue
                for model in coeffs:
                    try: coef,_,_,_,matrix=fit_surface(zb,vz,db,vd,xe,model)
                    except (ValueError,np.linalg.LinAlgError): continue
                    k=np.exp(np.clip(matrix@coef,-3,3))
                    for target,values in predictions(zb,db,k).items(): boot[model][target].append(values)
            study[f'{cat}/{channel}']={}
            for target,indices in TARGETS.items():
                observed=z[:,indices].sum(axis=1); vo=vz[:,indices].sum(axis=1)
                dt=d[:,indices].sum(axis=1); vdt=vd[:,indices].sum(axis=1)
                alpha=dt/dc
                vp=alpha*alpha*vc+c*c*(vdt/dc**2+dt*dt*vdc/dc**4)
                estimates={'nominal':(nominal[target],np.diag(vp))}
                for model,(coef,matrix) in coeffs.items():
                    samples=np.array(boot[model][target])
                    if len(samples)<50: continue
                    estimates[model]=(predictions(z,d,np.exp(matrix@coef))[target],np.cov(samples,rowvar=False))
                fig,(ax,ar)=plt.subplots(2,1,figsize=(9,7),sharex=True,gridspec_kw={'height_ratios':[3,1]},constrained_layout=True)
                width=np.diff(xe)
                ax.errorbar(xc,observed/width,np.sqrt(vo)/width,fmt='o',color='black',label='Held-out Data - TT/WW - VZ')
                for model,(pred,covpred) in estimates.items():
                    totalcov=covpred+np.diag(vo)
                    residual=observed-pred
                    score=float(residual@np.linalg.pinv(totalcov)@residual)
                    norm=float(observed.sum()/pred.sum())
                    row=dict(category=cat,channel=channel,training_channel=channel,target=target,model=model,observed=observed.sum(),predicted=pred.sum(),observed_over_predicted=norm,stat_chi2=score,bins=nz,bootstrap_success=len(boot[model][target]) if model!='nominal' else 0)
                    summary.append(row)
                    for i in range(nz):
                        predrows.append(dict(category=cat,channel=channel,target=target,model=model,xlow=xe[i],xhigh=xe[i+1],observed=observed[i],observed_stat=np.sqrt(vo[i]),nominal=nominal[target][i],corrected=pred[i],corrected_stat=np.sqrt(covpred[i,i]),correction_factor=pred[i]/nominal[target][i]))
                    study[f'{cat}/{channel}'].setdefault(target,{})[model]=row
                    line=ax.stairs(pred/width,xe,label=model.replace('_',' '))
                    color=line.get_edgecolor()
                    band=np.sqrt(np.diag(covpred))/np.maximum(pred,1e-20)
                    ar.plot(xc,observed/pred,'o-',color=color,label=model)
                    ar.fill_between(xc,1-band,1+band,color=color,alpha=.07)
                ar.axhline(1,color='gray',ls='--'); ar.set_ylim(.5,1.5)
                ar.set_xlabel('m(2l2q) [GeV]'); ar.set_ylabel('Observed / prediction')
                ax.set_ylabel('Events / GeV'); ax.legend(fontsize=8)
                ax.set_title(f'2018 {cat} untag, {channel}: held-out {target}\nFits use only 40-70, 105-135, 165-180 GeV hadronic-mass sidebands')
                savefig(fig,args.output,f'{cat}_{channel}_{target}_corrections')
        # Train on one lepton channel's sidebands, then validate on the other
        # channel. No target-region observations enter calibration.
        for source,test in [('ee','mumu'),('mumu','ee')]:
            zs,vzs,ds,vds=channels[source][:4]
            zt,vzt,dt,vdt=channels[test][:4]
            for model in ['separable','local_interpolation']:
                try: coef,_,_,_,matrix=fit_surface(zs,vzs,ds,vds,xe,model)
                except (ValueError,np.linalg.LinAlgError):continue
                central=predictions(zt,dt,np.exp(matrix@coef))
                samples={target:[] for target in TARGETS}
                for b in range(args.bootstrap):
                    zsb=rng.normal(zs,np.sqrt(vzs));dsb=rng.normal(ds,np.sqrt(vds))
                    ztb=rng.normal(zt,np.sqrt(vzt));dtb=rng.normal(dt,np.sqrt(vdt))
                    if np.any(dtb[:,TRAIN].sum(axis=1)<=0):continue
                    if any(np.any(dtb[:,ix].sum(axis=1)<=0) for ix in TARGETS.values()):continue
                    try:cb,_,_,_,mb=fit_surface(zsb,vzs,dsb,vds,xe,model)
                    except (ValueError,np.linalg.LinAlgError):continue
                    for target,values in predictions(ztb,dtb,np.exp(np.clip(mb@cb,-3,3))).items():samples[target].append(values)
                for target,ix in TARGETS.items():
                    if len(samples[target])<50:continue
                    obs=zt[:,ix].sum(axis=1);var=vzt[:,ix].sum(axis=1)
                    covariance=np.cov(np.array(samples[target]),rowvar=False)+np.diag(var)
                    residual=obs-central[target]
                    summary.append(dict(category=cat,channel=test,training_channel=source,target=target,model=model,observed=obs.sum(),predicted=central[target].sum(),observed_over_predicted=obs.sum()/central[target].sum(),stat_chi2=float(residual@np.linalg.pinv(covariance)@residual),bins=len(obs),bootstrap_success=len(samples[target])))
    for filename,rows in [('cells.csv',cells),('double_ratios.csv',drrows),('region_double_ratios.csv',regionrows),('double_ratio_tests.csv',drtests),('validation.csv',summary),('fit_quality.csv',fitrows),('correction_factors.csv',predrows)]:
        with (args.output/filename).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (args.output/'summary.json').write_text(json.dumps(study,indent=2))
    for row in summary:
        if row['channel']=='combined':print(row,flush=True)


if __name__=='__main__':main()
