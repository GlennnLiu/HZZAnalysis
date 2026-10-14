"""Read-only output checks for yields, correlations, overflow and plot coverage."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', type=Path)
    args = p.parse_args()
    base = args.directory
    rows = list(csv.DictReader((base/'yields.csv').open()))
    audit = json.loads((base/'input_audit.json').read_text())
    meta = json.loads((base/'metadata.json').read_text())
    assert len(rows) == 36
    audit_by_name = {Path(a['path']).name: a for a in audit}
    with np.load(base/'histograms.npz') as f:
        for r in rows:
            prefix = '__'.join([r['category'], r['tag'], r['channel'], r['target']])+'__'
            def get(name): return f[prefix+name]
            assert get('edges')[0] == 300 and np.isinf(get('edges')[-1])
            for name in ['pred_nozh', 'pred_withzh', 'pred_total_withzh', 'obs_nozh', 'obs_withzh']:
                value, variance = get(name)
                assert np.all(np.isfinite(value)) and np.all(variance >= 0)
                np.testing.assert_allclose(value.sum(), float(r[name]), rtol=1e-10, atol=1e-8)
                np.testing.assert_allclose(np.sqrt(variance.sum()), float(r[name+'_stat']), rtol=1e-10, atol=1e-8)
            np.testing.assert_allclose(get('pred_nozh')[0]-get('pred_withzh')[0], get('transferred_zh')[0], atol=1e-8)
            np.testing.assert_allclose(get('pred_nozh')[1]+get('pred_withzh')[1]-2*get('cov_before_after_prediction'),
                                       get('delta_dy_prediction')[1], atol=1e-8)
            if r['channel'] != 'combined':
                for reg in ['CR', r['target']]:
                    for proc in ['Data', 'DY', 'TTplusWW', 'VZ', 'ZHbb']:
                        name = '{}_{}_{}_{}_2018_{}.root'.format(proc,r['channel'],r['category'],r['tag'],reg)
                        key = ('CR' if reg == 'CR' else 'target')+'_'+proc
                        np.testing.assert_allclose(get(key)[0].sum(), audit_by_name[name]['selected_weight'],
                                                   rtol=1e-10, atol=1e-8)
            else:
                for key in ['pred_nozh','pred_withzh','obs_nozh','obs_withzh','transferred_zh']:
                    channel_sum = sum(f['__'.join([r['category'],r['tag'],ch,r['target']])+'__'+key]
                                      for ch in ['ee','mumu'])
                    np.testing.assert_allclose(get(key), channel_sum, rtol=1e-10, atol=1e-8)
    pngs = len(list(base.glob('*.png'))); pdfs = len(list(base.glob('*.pdf')))
    assert pngs+pdfs == meta['plots'], (pngs, pdfs, meta['plots'])
    report = dict(status='passed', target_rows=len(rows), cached_sources=len(audit),
                  pngs=pngs, pdfs=pdfs, all_search_mass_and_tail_preserved=True,
                  before_after_identities=True, correlated_difference_variance=True,
                  combined_channel_sum=True)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__': main()
