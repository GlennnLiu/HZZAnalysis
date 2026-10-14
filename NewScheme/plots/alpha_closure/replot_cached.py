"""Replot saved alpha closure predictions; do not reread trees or refit."""
from pathlib import Path
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from plot_helpers import draw,plt,hep,TAGS


if __name__=='__main__':
    root=Path(__file__).resolve().parents[3]/'Plot/alpha_closure_pnet_2018'; plt.style.use(hep.style.CMS)
    with np.load(root/'histograms.npz') as archive:
        for region in ('SR','VR'):
            with PdfPages(root/f'all_{region}_alpha_DataBkg_PNet.pdf') as bundle:
                for channel in ('ee','mumu'):
                    fig=plt.figure(figsize=(33,20)); outer=fig.add_gridspec(2,3,left=.065,right=.99,bottom=.09,top=.95,hspace=.40,wspace=.35)
                    for row,category in enumerate(('resolved','merged')):
                        for col,tag in enumerate(TAGS):
                            key=f'{category}_{tag}_{channel}_{region}'
                            edges=archive[key+'__edges']; pred=archive[key+'__prediction']; data=archive[key+'__data']
                            zero=np.zeros_like(pred)
                            parts=dict(DY=pred,TTplusWW=zero,VZ=zero,Data=np.tile(data[0],(5,1)))
                            kwargs=dict(edges_override=edges,channel=channel,region=region)
                            one=draw(parts,category,tag,'m2l2q',**kwargs)
                            for ext in ('png','pdf'): one.savefig(root/f'{key}_alpha_DataBkg_PNet.{ext}',dpi=140)
                            bundle.savefig(one); plt.close(one)
                            sub=outer[row,col].subgridspec(2,1,height_ratios=[3.3,1],hspace=.05)
                            ax=fig.add_subplot(sub[0]); lower=fig.add_subplot(sub[1],sharex=ax)
                            draw(parts,category,tag,'m2l2q',(ax,lower),**kwargs)
                    for ext in ('png','pdf'): fig.savefig(root/f'six_categories_{channel}_{region}_alpha_DataBkg_PNet.{ext}',dpi=120)
                    plt.close(fig)
            print('Replotted '+region,flush=True)
