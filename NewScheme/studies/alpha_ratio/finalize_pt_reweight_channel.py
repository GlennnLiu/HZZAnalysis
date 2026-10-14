"""Check the complete channel-dependent package and refresh its ZIP."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    args = parser.parse_args()
    folder = args.study.resolve()
    metadata = json.loads((folder/'metadata.json').read_text())
    assert len(list(folder.glob('*.png'))) == 147
    assert len(list(folder.glob('*.pdf'))) == 150
    assert all(path.with_suffix('.pdf').is_file() for path in folder.glob('*.png'))
    assert metadata['bundle_pages'] == {'mass_and_factors': 93, 'SR_pt': 18, 'CR_pt': 36}
    assert len(metadata['fits']) == 24
    assert len(metadata['calibration_checks']) == 36
    assert len(metadata['mass_pt_projection_checks']) == 36
    assert not metadata['shared_tag_corrections'] and not metadata['ZH_subtracted']
    digest = hashlib.sha256((folder/'source_histograms.npz').read_bytes()).hexdigest()
    assert digest == metadata['source_sha256']
    with np.load(folder/'results.npz') as data:
        for key in data.files:
            if '__calibration__' in key and key.endswith('__dy_reweighted'):
                root = key[:-len('dy_reweighted')]
                np.testing.assert_allclose(data[key][0], data[root+'observed'][0], rtol=1e-12, atol=1e-9)
                np.testing.assert_allclose(data[key+'_closure_bins'][0], 1., rtol=0, atol=1e-12)
                np.testing.assert_allclose(data[key+'_closure_bins'][1], 0., rtol=0, atol=1e-12)
    result = subprocess.run([sys.executable, '-m', 'unittest', 'test_pt_reweight',
                             'test_pt_reweight_channel', '-v'], cwd=str(folder),
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
    (folder/'tests.log').write_text(result.stdout)
    document = Path(__file__).parent/'PT_CHANNEL_STUDY.md'
    if document.exists():
        shutil.copy2(document, folder/'README.md')
    if Path(__file__).resolve().parent != folder:
        shutil.copy2(Path(__file__), folder/Path(__file__).name)
    validation = dict(status='passed', regression_tests=13, calibration_spectra=36,
                      pt_mass_projections=36, independent_channel_fits=24,
                      individual_png_pdf_pairs=147, PDF_files=150, source_sha256=digest,
                      additionally_merged_mass_grids=metadata['additionally_merged_mass_grids'])
    (folder/'validation.json').write_text(json.dumps(validation, indent=2))
    metadata['script_sha256_all'] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                     for path in folder.glob('*.py')}
    (folder/'metadata.json').write_text(json.dumps(metadata, indent=2))
    target = folder.with_suffix('.zip')
    with zipfile.ZipFile(target, 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(folder.iterdir()):
            if path.is_file():
                archive.write(path, arcname=folder.name+'/'+path.name)
    print(json.dumps(dict(validation, ZIP_bytes=target.stat().st_size)), flush=True)


if __name__ == '__main__':
    main()
