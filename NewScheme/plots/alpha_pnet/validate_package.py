"""Read-only checks of saved alpha comparisons."""
import csv
from pathlib import Path
import numpy as np


if __name__=='__main__':
    root=Path(__file__).resolve().parents[3]/'Plot/alpha_pnet_comparison_2018'
    with np.load(root/'histograms.npz') as archive:
        assert len(archive.files)==192
        for region in ('SR','VR'):
            for channel in ('ee','mumu'):
                for category in ('resolved','merged'):
                    for tag in ('untag','btag','vbftag'):
                        key=f'{category}_{tag}_{channel}_{region}'
                        e=archive[key+'__edges']; t=archive[key+'__target']; c=archive[key+'__source']
                        v=archive[key+'__alpha']; valid=archive[key+'__valid']
                        assert e[0]==300 and np.all(np.diff(e)>0)
                        assert t.shape==c.shape==(5,len(e)-1)
                        assert v.shape==(8,len(e)-1)
                        assert np.isfinite(v[:,valid]).all()
                        np.testing.assert_allclose(v[0,valid],t[0,valid]/c[0,valid])
                        np.testing.assert_allclose(v[1,valid],t[1,valid]/c[1,valid])
                        np.testing.assert_allclose(v[5,valid],v[1,valid]/v[0,valid])
                        assert np.all(v[2:4,valid]>=0) and np.all(v[6:,valid]>=0)
                        assert np.all(v[4,valid]**2<=v[2,valid]*v[3,valid]+1e-10)
                        if category=='resolved':
                            np.testing.assert_array_equal(v[0,valid],v[1,valid])
                            np.testing.assert_allclose(v[6:,valid],0,atol=1e-12)
                        for ext in ('png','pdf'): assert (root/f'{key}_alpha_PNet.{ext}').stat().st_size>0
                for ext in ('png','pdf'): assert (root/f'six_categories_{channel}_{region}_alpha_PNet.{ext}').stat().st_size>0
            assert (root/f'all_{region}_alpha_PNet.pdf').stat().st_size>0
    assert len(list(root.glob('*.png')))==28
    assert len(list(root.glob('*.pdf')))==30
    with (root/'summary.csv').open() as handle: assert len(list(csv.DictReader(handle)))==36
    print('PASS: 24 paired-alpha comparisons, original edges, nonnegative variances, exact resolved cancellation, 28 PNGs / 30 PDFs')
