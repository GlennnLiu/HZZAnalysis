"""Readable nominal-versus-candidate figures from the double-ratio tables."""
import argparse
import csv
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    args=p.parse_args()
    with (args.directory/'correction_factors.csv').open() as f:rows=list(csv.DictReader(f))
    for cat in ['resolved','merged']:
        for target in ['SR','VR']:
            selected=[r for r in rows if r['category']==cat and r['target']==target and r['channel']=='combined']
            records={model:[r for r in selected if r['model']==model] for model in ['nominal','separable']}
            if not all(records.values()):continue
            rr=records['nominal'];lo=np.array([float(r['xlow']) for r in rr]);hi=np.array([float(r['xhigh']) for r in rr])
            x=(lo+hi)/2;edges=np.r_[lo,hi[-1]];width=hi-lo
            obs=np.array([float(r['observed']) for r in rr]);error=np.array([float(r['observed_stat']) for r in rr])
            fig,(ax,ar)=plt.subplots(2,1,figsize=(8,7),sharex=True,gridspec_kw={'height_ratios':[3,1]},constrained_layout=True)
            ax.errorbar(x,obs/width,error/width,fmt='o',color='black',label='Data - TT/WW - VZ')
            for model,color,label,marker in [('nominal','tab:blue','Nominal alpha','o'),('separable','tab:orange','Sideband correction candidate','^')]:
                pred=np.array([float(r['corrected']) for r in records[model]])
                err=np.array([float(r['corrected_stat']) for r in records[model]])
                ax.stairs(pred/width,edges,color=color,label=label)
                ar.errorbar(x,obs/pred,error/pred,fmt=marker,color=color,label=label)
                band=err/pred
                ar.fill_between(edges,np.r_[1-band,(1-band)[-1]],np.r_[1+band,(1+band)[-1]],step='post',color=color,alpha=.12)
            ax.set_yscale('log');ax.set_ylabel('Events / GeV');ax.legend(fontsize=9)
            ar.set_xscale('log');ar.set_xticks([300,400,600,1000,1600,3000]);ar.set_xticklabels(['300','400','600','1000','1600','3000'])
            ar.axhline(1,color='gray',ls='--');ar.set_ylim(.6,1.4)
            ar.set_ylabel('Data / prediction');ar.set_xlabel('m(2l2q) [GeV]')
            ax.set_title(f'2018 {cat} untag, ee + mumu: {target}\nCorrection fitted only in hadronic-mass sidebands')
            for ext in ['png','pdf']:fig.savefig(args.directory/f'{cat}_combined_{target}_candidate_summary.{ext}',dpi=180)
            plt.close(fig)


if __name__=='__main__':main()
