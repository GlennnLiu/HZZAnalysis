"""Independent final audit of all signal augmentation and reduction products."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[key]='1'
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import uproot
import produce_signal_reduced as producer

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--augmented',type=Path,required=True); p.add_argument('--reduced',type=Path,required=True)
    args=p.parse_args(); here=Path(__file__).resolve().parent
    manifest=json.loads((args.augmented/'augmentation_manifest.json').read_text())
    assert len(manifest)==24 and all(r['status']=='passed' for r in manifest)
    expected_samples={(era,f'{mode}{mass}') for era in ('2016preAPV','2016postAPV','2017','2018') for mode in ('ggh','vbf') for mass in (500,1000,3000)}
    assert {(r['era'],r['sample']) for r in manifest}==expected_samples
    for r in manifest:
        src=Path(r['source']); assert (src.stat().st_size,src.stat().st_mtime_ns)==(r['source_size'],r['source_mtime_ns'])
        assert hashlib.sha256(Path(r['calibration']).read_bytes()).hexdigest()==r['calibration_sha256']
        with uproot.open(r['output']) as f:
            t=f['passedEvents']; assert t.num_entries==r['entries']
            a=t.arrays(['w_deepjet','w_deepjet_up','w_deepjet_dn','w_deepjet_valid'],library='np'); valid=a['w_deepjet_valid']
            assert int(valid.sum())==r['selected']
            for field in ('w_deepjet','w_deepjet_up','w_deepjet_dn'):
                assert np.isfinite(a[field]).all() and (a[field][~valid]==-999).all()
            assert (a['w_deepjet'][valid]>=0).all()
    files=[]; sum_entries=0; corrections={}
    for era in ('2016preAPV','2016postAPV','2017','2018'):
        for cat in ('resolved','merged'):
            label=f'{era}_{cat}'; path=args.reduced/f'validation_{label}.json'
            assert (args.reduced/f'{label}.done').is_file()
            report=json.loads(path.read_text()); assert report['status']=='passed' and not report['pilot'] and len(report['outputs'])==108
            cp=Path(report['config']); assert hashlib.sha256(cp.read_bytes()).hexdigest()==report['config_sha256']
            original=here.parent/'reduced_20261009_pnet'/'cards'/cp.name
            assert hashlib.sha256(original.read_bytes()).hexdigest()==report['config_sha256']
            for row in report['outputs']:
                file=args.reduced/row['file']; files.append(file.name); totals=np.zeros(3); n=0
                with uproot.open(file) as f:
                    assert set(f['tree'].keys())==set(report['branches'])
                    for a in f['tree'].iterate(step_size=100000,library='np'):
                        producer.validate(a,cat,row['tag'],row['region'],int(row['sample'][3:])); n+=len(a['weight'])
                        totals+=np.array([a['weight'].sum(),(a['weight']**2).sum(),a['weight_no_pnet'].sum()])
                assert n==row['entries']; sum_entries+=n
                np.testing.assert_allclose(totals,[row[k] for k in ('sumw','sumw2','sumw_no_pnet')],rtol=1e-10,atol=1e-8)
            corrections[label]=dict(files=len(report['outputs']),entries=sum(r['entries'] for r in report['outputs']))
    assert len(files)==864 and len(set(files))==864 and len(list(args.reduced.glob('*.root')))==864
    assert not list(args.reduced.glob('*.partial'))
    result=dict(status='passed',augmented_files=24,reduced_files=864,stored_entries_including_VR_CR_overlap=sum_entries,era_topology=corrections,
                original_ntuples_unchanged=True,frozen_cut_yamls_unchanged=True,normalization='1 pb benchmark per mass point; original era luminosities')
    (args.reduced/'complete_validation.json').write_text(json.dumps(result,indent=2)); print(json.dumps(result,indent=2))
if __name__=='__main__': main()
