"""Read-only checks of the saved full-range score histograms and selection audit."""
from pathlib import Path
import csv
import json
import hashlib
import numpy as np

p=Path(__file__).resolve().parents[3]/'Plot/run2_PNet_fullrange'
meta=json.loads((p/'metadata.json').read_text())
assert meta['SR_Data_used'] is False and meta['luminosity_display']==138
assert len(meta['source_files'])==44 and len(meta['config_snapshots'])==4
for snapshot in meta['config_snapshots']:
    path=Path(snapshot['path'])
    if path.is_file():
        assert hashlib.sha256(path.read_bytes()).hexdigest()==snapshot['sha256'], 'Input YAML changed'
rows=list(csv.DictReader((p/'input_yields.csv').open()))
assert len(rows)==84 and not any(r['process']=='Data' and r['region']=='SR' for r in rows)
with np.load(p/'histograms.npz') as a:
    assert not any('_SR_' in k and '__Data' in k for k in a.files)
    for region in ('SR','CR'):
        prefix=f'merged_{region}_particleNetZvsQCD'
        np.testing.assert_allclose(a[prefix+'__edges'],np.linspace(0,1,21))
        for proc in ('VZ','TTplusWW','DY')+(() if region=='SR' else ('Data',)):
            subset=[r for r in rows if r['region']==region and r['process']==proc]
            assert len(subset)=={'VZ':8,'TTplusWW':12,'DY':20,'Data':4}[proc]
            for r in subset:
                assert int(r['selected'])==int(r['valid'])+int(r['invalid'])
                assert int(r['valid'])==int(r['below_wp'])+int(r['tagged'])
            assert sum(int(r['below_wp']) for r in subset)>0
            for suffix,fieldprefix in (('',''),('__tagged_check','tagged_')):
                h=a[prefix+'__'+proc+suffix]
                assert np.isfinite(h).all() and (h[1]>=0).all()
                np.testing.assert_allclose(h.sum(axis=1),[sum(float(r[fieldprefix+k]) for r in subset) for k in ('sumw','sumw2')],rtol=1e-10,atol=1e-6)
                if proc=='Data': np.testing.assert_array_equal(h[0],h[1])
            if region=='CR':
                assert sum(int(r['bridge_entries']) for r in subset)>0
                assert sum(int(r['VR_entries']) for r in subset)>0
        for ext in ('png','pdf'): assert (p/f'merged_{region}_particleNetZvsQCD_fullrange.{ext}').is_file()
for ext in ('png','pdf'): assert (p/f'SR_CR_fullrange_particleNetZvsQCD.{ext}').is_file()
assert len(list(p.glob('*.png')))==3 and len(list(p.glob('*.pdf')))==4
cr=[r for r in rows if r['region']=='CR' and r['process']=='Data']
print('CR Data: full range',sum(int(r['valid']) for r in cr),'score<=0.9',sum(int(r['below_wp']) for r in cr),'score>0.9',sum(int(r['tagged']) for r in cr))
print('PASS: 44 unique input files, 84 sample/region audits, full-range yields and variances, bridge/VR retained, SR Data absent, all figures')
