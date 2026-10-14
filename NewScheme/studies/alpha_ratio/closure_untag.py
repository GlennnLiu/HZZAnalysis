"""2018 untagged alpha closure using internally consistent reduced trees."""
import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2


def read(base, process, channel, cat, region):
    with uproot.open(base / f'{process}_{channel}_{cat}_untag_2018_{region}.root') as f:
        return f['tree'].arrays(['mass', 'mass2q', 'weight'], library='np')


def hist(a, edges, region, window):
    q = a['mass2q']
    mask = np.ones(len(q), dtype=bool)
    if region == 'CR':
        mask &= ~((q > window[0]) & (q < window[1]))
    m, w = a['mass'][mask], a['weight'][mask]
    return np.histogram(m, edges, weights=w)[0], np.histogram(m, edges, weights=w*w)[0]


def ratio(n, d):
    return np.divide(n, d, out=np.full_like(n, np.nan), where=d > 0)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base', type=Path, required=True)
    p.add_argument('--cr-base', type=Path)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    edges = np.array([300,400,500,600,750,900,1100,1400,1800,2300,3500.])
    rows, bins = [], []
    audit = []
    for cat in ['resolved', 'merged']:
        for channel in ['ee', 'mumu']:
            arrays = {(proc, reg): read((args.cr_base or args.base) if reg=='CR' else args.base,proc,channel,cat,reg)
                      for proc in ['Data','DY','TTplusWW','VZ'] for reg in ['CR','SR','VR']}
            if args.cr_base:
                qe=np.array([40,55,70,85,105,120,135,150,165,180.])
                qh={}
                for proc in ['Data','DY','TTplusWW','VZ']:
                    qh[proc]=[np.zeros(len(qe)-1),np.zeros(len(qe)-1)]
                    for reg in ['CR','SR']:
                        a=arrays[proc,reg]
                        valid=(a['mass']>=300)&(a['mass']<3500)
                        q,w=a['mass2q'][valid],a['weight'][valid]
                        qh[proc][0]+=np.histogram(q,qe,weights=w)[0]
                        qh[proc][1]+=np.histogram(q,qe,weights=w*w)[0]
                z=qh['Data'][0]-qh['TTplusWW'][0]-qh['VZ'][0]
                vz=qh['Data'][1]+qh['TTplusWW'][1]+qh['VZ'][1]
                d,vd=qh['DY']
                x=(qe[:-1]+qe[1:])/2
                fig,(ax,ar)=plt.subplots(2,1,figsize=(8,7),sharex=True,gridspec_kw={'height_ratios':[3,1]},constrained_layout=True)
                ax.errorbar(x,z/np.diff(qe),np.sqrt(vz)/np.diff(qe),fmt='o',color='black',label='Data − TT/WW − VZ')
                ax.stairs(d/np.diff(qe),qe,label='DY MC',color='tab:orange')
                ax.set_ylabel('Events / GeV'); ax.legend()
                ax.set_title(f'2018 {cat}, {channel}, untag; 300 < m(2ℓ2q) < 3500 GeV')
                r=ratio(z,d)
                vr=ratio(vz,d*d)+ratio(z*z*vd,d**4)
                ar.errorbar(x,r,np.sqrt(vr),fmt='o',color='black')
                ar.axhline(1,color='gray',ls='--'); ar.set_ylim(.5,1.6)
                ar.set_ylabel('Signed data / DY'); ar.set_xlabel('m(qq) [GeV]')
                for a in [ax,ar]:
                    a.axvspan(70,105,color='tab:green',alpha=.1)
                    a.axvspan(135,165,color='tab:purple',alpha=.1)
                for ext in ['png','pdf']: fig.savefig(args.output/f'{cat}_{channel}_mass2q.{ext}',dpi=180)
                plt.close(fig)
            for reg in ['CR','SR','VR']:
                q = arrays['Data',reg]['mass2q']
                audit.append(f'{cat}/{channel}/{reg}: entries={len(q)}, mass2q=[{q.min():g},{q.max():g}], SR_fraction={np.mean((q>70)&(q<105)):.6f}, H_fraction={np.mean((q>105)&(q<135)):.6f}, VR_fraction={np.mean((q>135)&(q<165)):.6f}')
            for target in ['SR','VR']:
                for window in [(135,165),(105,165)]:
                    hs = {key: hist(a,edges,key[1],window) for key,a in arrays.items()}
                    def signed(reg):
                        return (hs['Data',reg][0]-hs['TTplusWW',reg][0]-hs['VZ',reg][0],
                                hs['Data',reg][1]+hs['TTplusWW',reg][1]+hs['VZ',reg][1])
                    c, vc = signed('CR')
                    t, vt = signed(target)
                    dc, vdc = hs['DY','CR']
                    dt, vdt = hs['DY',target]
                    alpha = ratio(dt, dc)
                    va = ratio(vdt, dc*dc)+ratio(dt*dt*vdc,dc**4)
                    pred = alpha*c
                    vp = alpha*alpha*vc+c*c*va
                    valid = np.isfinite(pred)&np.isfinite(vp)&((vt+vp)>0)&(dc>0)
                    x2 = np.sum((t[valid]-pred[valid])**2/(vt[valid]+vp[valid]))
                    n = int(valid.sum())
                    rawx2 = np.sum((t-dt)**2/(vt+vdt))
                    row = dict(category=cat,channel=channel,target=target,exclude=f'{window[0]}-{window[1]}',observed=t.sum(),observed_stat=np.sqrt(vt.sum()),predicted=np.nansum(pred),prediction_stat=np.sqrt(np.nansum(vp)),dy=dt.sum(),obs_over_pred=t.sum()/np.nansum(pred),chi2=x2,ndof=n,pvalue=chi2.sf(x2,n),raw_dy_chi2=rawx2)
                    rows.append(row)
                    for i in range(len(t)):
                        bins.append(dict(category=cat,channel=channel,target=target,exclude=row['exclude'],low=edges[i],high=edges[i+1],observed=t[i],observed_variance=vt[i],predicted=pred[i],prediction_variance=vp[i],dy=dt[i],dy_variance=vdt[i],cr=c[i],dy_cr=dc[i],alpha=alpha[i],closure=ratio(t,pred)[i]))
                    fig,(ax,ar)=plt.subplots(2,1,figsize=(8,7),sharex=True,gridspec_kw={'height_ratios':[3,1]},constrained_layout=True)
                    x=(edges[:-1]+edges[1:])/2
                    widths=np.diff(edges)
                    ax.errorbar(x,t/widths,np.sqrt(vt)/widths,fmt='o',color='black',label='Data − TT/WW − VZ')
                    ax.stairs(pred/widths,edges,color='tab:blue',label='Nominal DY alpha × signed CR')
                    ax.fill_between(x,(pred-np.sqrt(vp))/widths,(pred+np.sqrt(vp))/widths,color='tab:blue',alpha=.2,label='Prediction statistical uncertainty')
                    ax.stairs(dt/widths,edges,color='tab:orange',label='Direct DY MC')
                    ax.set_ylabel('Events / GeV')
                    ax.set_title(f'2018 {cat}, {channel}, untag — {target}\nCR excludes {window[0]}–{window[1]} GeV; 300 < mass < 3500 GeV')
                    ax.legend(fontsize=9)
                    r=ratio(t,pred)
                    ar.errorbar(x,r,ratio(np.sqrt(vt),pred),fmt='o',color='black')
                    band=ratio(np.sqrt(vp),pred)
                    ar.fill_between(x,1-band,1+band,color='tab:blue',alpha=.2)
                    ar.axhline(1,color='gray',ls='--')
                    ar.set_ylim(0,2)
                    ar.set_ylabel('Data / alpha')
                    ar.set_xlabel('m(2ℓ2q) [GeV]')
                    stem=f'{cat}_{channel}_untag_{target}_exclude_{window[0]}_{window[1]}'
                    for ext in ['png','pdf']:
                        fig.savefig(args.output/f'{stem}.{ext}',dpi=180)
                    plt.close(fig)
    for filename, records in [('yields.csv',rows),('bins.csv',bins)]:
        with (args.output/filename).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
    (args.output/'input_audit.txt').write_text('\n'.join(audit)+'\n')
    for line in audit: print(line)
    for r in rows: print(r)


if __name__ == '__main__':
    main()
