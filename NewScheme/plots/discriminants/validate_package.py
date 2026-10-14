"""Validate the retained Run-2 discriminant plot package, without ROOT reads."""
import csv
import json
from pathlib import Path
import numpy as np

p=Path(__file__).resolve().parents[3]/'Plot/run2_discriminants_PNet'
meta=json.loads((p/'metadata.json').read_text())
assert meta['SR_Data_read'] is False and meta['luminosity_display']==138
assert len(meta['source_files'])==336 and len(meta['validation_report_sha256'])==8
assert not any('/Data_' in s['file'] and '_SR.root' in s['file'] for s in meta['source_files'])
rows=list(csv.DictReader((p/'validity_audit.csv').open()))
assert len(rows)==1008
with np.load(p/'histograms.npz') as a:
    assert not any('_SR_' in k and '__Data' in k for k in a.files)
    for cat in ('resolved','merged'):
        for region in ('SR','CR'):
            for var in ('kd','kd_jjvbf','kd_jvbf'):
                key=f'{cat}_{region}_{var}'
                np.testing.assert_allclose(a[key+'__edges'],np.linspace(0,1,21))
                for proc in ('VZ','TTplusWW','DY')+(() if region=='SR' else ('Data',)):
                    matching=[r for r in rows if (r['category'],r['region'],r['variable'],r['process'])==(cat,region,var,proc)]
                    assert len(matching)==24
                    for r in matching:
                        assert int(r['selected'])==sum(int(r[k]) for k in ('valid','negative','nonfinite','above_one'))
                    h=a[key+'__'+proc]
                    assert np.isfinite(h).all() and (h[1]>=0).all()
                    np.testing.assert_allclose(h.sum(axis=1),[sum(float(r[k]) for r in matching) for k in ('sumw','sumw2')],rtol=1e-10,atol=1e-6)
                    if proc=='Data': np.testing.assert_array_equal(h[0],h[1])
                for ext in ('png','pdf'): assert (p/f'{key}.{ext}').is_file()
            for ext in ('png','pdf'): assert (p/f'{cat}_{region}_three_discriminants.{ext}').is_file()
assert len(list(p.glob('*.png')))==16 and len(list(p.glob('*.pdf')))==17
summary={}
for cat in ('resolved','merged'):
    for var in ('kd','kd_jjvbf','kd_jvbf'):
        subset=[r for r in rows if r['category']==cat and r['region']=='CR' and r['process']=='Data' and r['variable']==var]
        summary[f'{cat}_CR_Data_{var}']={k:sum(int(r[k]) for r in subset) for k in ('selected','valid','negative','nonfinite','above_one','both_vbf_valid')}
print(json.dumps(summary,indent=2))
print('PASS: all valid yields and variances, 336 input files, SR Data absent, 12 distributions and 4 overview pages')
