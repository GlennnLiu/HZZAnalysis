"""Read-only numerical and completeness checks after rendering."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='1'
import json
from pathlib import Path
import numpy as np
from prepare_histograms import MC,TAGS,CHANNELS,SIGNALS
from render_plots import mass_parts,read_npz

base=Path('/cms/user/liugeliang/CMS/HZZAnalysis/Plot')
cache=base/'signal_overlay_Run2_20261013'
metadata=json.loads((cache/'metadata.json').read_text())
assert metadata['SR_Data_read'] is False
assert metadata['signal_xsec_pb']==1
mass=read_npz(cache/'mass_histograms.npz')
sig=read_npz(cache/'signal_histograms.npz')
full=read_npz(cache/'signal_fullrange_histograms.npz')
for cat in ('resolved','merged'):
    for p in MC+SIGNALS:
        inclusive=mass_parts(mass,cat,'all','combined')[p]
        expected=sum((mass[f'{cat}__{t}__{c}__{p}'] for t in TAGS for c in CHANNELS),np.zeros_like(inclusive))
        np.testing.assert_allclose(inclusive,expected)
for region in ('SR','CR'):
    for p in SIGNALS:
        key=f'merged_{region}_particleNetZvsQCD__{p}'
        np.testing.assert_allclose(full[key+'__tagged_check'],sig[key+'_no_pnet'],rtol=1e-10,atol=1e-7)
for folder in ('run2_discriminants_PNet','run2_PNet_scores','run2_PNet_fullrange'):
    bg=read_npz(base/folder/'histograms.npz')
    assert not any('_SR_' in k and '__Data' in k for k in bg)
    destination=base/(folder+'_signals')
    assert {p.name for p in (base/folder).glob('*.pdf')}=={p.name for p in destination.glob('*.pdf')}
expected={'run2_discriminants_PNet_signals':(16,17),'run2_PNet_scores_signals':(6,7),
          'run2_PNet_scores_linear_signals':(6,7),'run2_PNet_fullrange_signals':(3,4),
          'run2_PNet_fullrange_linear_signals':(3,4),'run2_SR_m2l2q_signals':(28,29)}
for name,counts in expected.items():
    folder=base/name
    assert (len(list(folder.glob('*.png'))),len(list(folder.glob('*.pdf'))))==counts,(name,counts)
    assert all(p.stat().st_size>1000 for ext in ('*.png','*.pdf') for p in folder.glob(ext))
print('PASS: 62 PNGs, 68 PDFs; channel/tag sums and weighted variances; full-range regression; original plot names retained; SR Data absent.')
