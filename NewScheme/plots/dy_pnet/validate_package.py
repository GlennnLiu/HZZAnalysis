"""Read-only final checks and concise yield summary."""
import csv
from pathlib import Path
import numpy as np
from plot_dy_pnet import ERAS,TAGS,ratio


if __name__=='__main__':
    root=Path(__file__).resolve().parents[3]/'Plot/dy_pnet_comparison_allEras'
    archive=np.load(root/'histograms.npz')
    assert len(archive.files)==192
    for era in ERAS:
        for category in ('resolved','merged'):
            for tag in TAGS:
                for kind in ('mqq','CR','SR','VR'):
                    key=f'{era}_{category}_{tag}_{kind}'
                    h=archive[key+'__hist']; edges=archive[key+'__edges']
                    assert h.shape==(5,len(edges)-1) and np.isfinite(h).all()
                    assert (h[2:4]>=0).all()
                    ratio(h)
                    if category=='resolved':
                        np.testing.assert_array_equal(h[0],h[1])
                        np.testing.assert_array_equal(h[2],h[3])
                        np.testing.assert_array_equal(h[2],h[4])
                    for ext in ('png','pdf'): assert (root/f'{key}.{ext}').stat().st_size>0
        for kind in ('mqq','CR','SR','VR'):
            for ext in ('png','pdf'): assert (root/f'six_categories_{era}_{kind}_DY_ParticleNet.{ext}').stat().st_size>0
        assert (root/f'all_{era}_DY_ParticleNet.pdf').stat().st_size>0
    assert len(list(root.glob('*.png')))==112
    assert len(list(root.glob('*.pdf')))==116
    with (root/'yields.csv').open() as handle: rows=list(csv.DictReader(handle))
    assert len(rows)==96
    for row in rows:
        if row['era']=='2018' and row['category']=='merged' and row['distribution']=='SR':
            print(f"2018 merged {row['tag']} SR: SF/noSF={float(row['ratio']):.6f}, yield shift={100*(float(row['ratio'])-1):+.2f}%")
    print('PASS: 96 individual comparisons, 16 overview pages, 4 PDF bundles; resolved histograms exactly unchanged')
