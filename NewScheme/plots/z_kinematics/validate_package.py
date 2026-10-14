"""Validate sums across eras/tags/channels, coverage, blinding and figures."""
import csv
import json
from pathlib import Path
import numpy as np


if __name__=='__main__':
    root=Path(__file__).resolve().parents[3]/'Plot/run2_z_kinematics_PNet'
    meta=json.loads((root/'metadata.json').read_text())
    assert meta['SR_Data_read'] is False and meta['luminosity_display']==138
    assert len(meta['source_files'])==336 and len(meta['validation_report_sha256'])==8
    assert not any('/Data_' in row['file'] and row['file'].endswith('_SR.root') for row in meta['source_files'])
    with (root/'input_yields.csv').open() as handle: rows=list(csv.DictReader(handle))
    assert len(rows)==336
    with np.load(root/'histograms.npz') as archive:
        for category in ('resolved','merged'):
            for region in ('SR','CR'):
                processes=('VZ','TTplusWW','DY') if region=='SR' else ('VZ','TTplusWW','DY','Data')
                for process in processes:
                    selected=[r for r in rows if r['category']==category and r['region']==region and r['process']==process]
                    assert len(selected)==24
                    expected=np.array([sum(float(r[k]) for r in selected) for k in ('sumw','sumw2')])
                    if region=='CR':
                        assert sum(int(r['bridge_entries']) for r in selected)>0
                        assert sum(int(r['VR_entries']) for r in selected)>0
                    for variable in ('pt2l','mass2l','pt2q','mass2q'):
                        key=f'{category}_{region}_{variable}'
                        h=archive[key+'__'+process]
                        assert np.isfinite(h).all() and np.all(h[1]>=0)
                        np.testing.assert_allclose(h.sum(axis=1),expected,rtol=1e-10,atol=1e-6)
                        if process=='Data': np.testing.assert_array_equal(h[0],h[1])
                        for ext in ('png','pdf'): assert (root/f'{key}.{ext}').stat().st_size>0
                if region=='SR': assert not any(k.startswith(category+'_SR_') and k.endswith('__Data') for k in archive.files)
                for ext in ('png','pdf'): assert (root/f'{category}_{region}_four_variables.{ext}').stat().st_size>0
    assert len(list(root.glob('*.png')))==20 and len(list(root.glob('*.pdf')))==21
    print('PASS: 16 kinematic distributions + 4 overview pages; exact full-input sums for all four variables; SR Data absent; CR bridge and VR retained')
