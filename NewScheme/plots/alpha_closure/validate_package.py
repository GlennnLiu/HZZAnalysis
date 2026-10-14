"""Check predictions, uncertainty propagation and figures; write a readable summary."""
import csv
from pathlib import Path
import numpy as np


if __name__=='__main__':
    root=Path(__file__).resolve().parents[3]/'Plot/alpha_closure_pnet_2018'; maximum=0.
    with np.load(root/'histograms.npz') as archive:
        for region in ('SR','VR'):
            for channel in ('ee','mumu'):
                for category in ('resolved','merged'):
                    for tag in ('untag','btag','vbftag'):
                        key=f'{category}_{tag}_{channel}_{region}'
                        a=archive[key+'__alpha']; z=archive[key+'__z_source']; o=archive[key+'__other_target']
                        d=archive[key+'__data']; p=archive[key+'__prediction']; mask=archive[key+'__chi2_mask']
                        valid=np.isfinite(p).all(axis=0)
                        np.testing.assert_array_equal(d[0],d[1])
                        np.testing.assert_allclose(p[:2,valid],a[:2,valid]*z[0,valid]+o[0,valid])
                        np.testing.assert_allclose(p[2:4,valid],a[:2,valid]**2*z[1,valid]+z[0,valid]**2*a[2:4,valid]+o[1,valid])
                        assert np.all(p[2,valid]+p[3,valid]-2*p[4,valid]>=-1e-6)
                        assert np.all(d[0,mask]>=5)
                        if category=='resolved': np.testing.assert_array_equal(p[0,valid],p[1,valid])
                        if mask.any(): maximum=max(maximum,float(np.max(np.abs(p[1,mask]-p[0,mask])/np.sqrt(p[2,mask]))))
                        for ext in ('png','pdf'): assert (root/f'{key}_alpha_DataBkg_PNet.{ext}').stat().st_size>0
    assert len(list(root.glob('*.png')))==28 and len(list(root.glob('*.pdf')))==30
    with (root/'compatibility.csv').open() as handle: rows=list(csv.DictReader(handle))
    assert len(rows)==72
    lines=['# Alpha-method SR/VR compatibility, 2018','',
           'Combined tests sum separate ee/mumu chi2 values; original channel binning. Statistical only.',
           '','| Region | Category | chi2/ndof no PNet | with PNet | Data/Bkg no PNet | with PNet |',
           '| --- | --- | ---: | ---: | ---: | ---: |']
    for region in ('SR','VR'):
        for category in ('resolved','merged'):
            for tag in ('untag','btag','vbftag'):
                selected=[r for r in rows if r['channel']=='combined_channel_tests' and r['region']==region and r['category']==category and r['tag']==tag]
                r0,r1=selected
                fmt=lambda r: f"{float(r['chi2']):.2f}/{r['ndof']}" if int(r['ndof']) else 'N/A (sparse)'
                lines.append(f"| {region} | {category} {tag} | {fmt(r0)} | {fmt(r1)} | {float(r0['data_over_bkg']):.4f} | {float(r1['data_over_bkg']):.4f} |")
    lines+=['','Yield ratios include only covered bins with defined alpha. Undefined-alpha Data yields are explicit in CSV;',
            'in merged VR they exclude 4 untag, 8 btag and 1 VBF events. No Gaussian chi2 is available for',
            'merged VR btag/VBF at original binning: CSV chi2=0, ndof=0 is an empty sum, not perfect closure.',
            '',f'Max SF shift / baseline prediction statistical sigma on chi2 bins: {maximum:.4f}.',
            '', 'No calibration/modeling systematics; alternative-method chi2 differences are not independent significance tests.']
    (root/'compatibility.md').write_text('\n'.join(lines)+'\n')
    print(f'PASS: 24 closure comparisons, propagated covariance, 28 PNGs / 30 PDFs; max shift/stat_sigma={maximum:.6f}')
