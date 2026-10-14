"""Validate plot completeness and package an existing pT-reweighting study."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--study',type=Path,required=True)
    args=parser.parse_args()
    figures=list(args.study.glob('*.png'))
    if len(figures)!=89:raise ValueError(f'Expected 89 individual figures, found {len(figures)}')
    for figure in figures:
        if not figure.with_suffix('.pdf').is_file():raise ValueError('Missing PDF for '+str(figure))
    for name in ['all_pt_reweight_plots.pdf','all_SR_pt_reweight.pdf','validation.json','README.md']:
        if not (args.study/name).is_file():raise ValueError('Missing '+name)
    for name in ('pt_reweight_core.py','test_pt_reweight.py','validate_pt_reweight.py','package_pt_reweight.py'):
        shutil.copy2(Path(__file__).parent/name,args.study/name)
    shutil.copy2(args.source/'source_histograms.npz',args.study/'source_histograms.npz')
    metadata=json.loads((args.study/'metadata.json').read_text())
    digest=hashlib.sha256((args.study/'source_histograms.npz').read_bytes()).hexdigest()
    if digest!=metadata['source_sha256']:raise ValueError('Source cache digest mismatch')
    metadata['script_sha256_all']={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in args.study.glob('*.py')}
    metadata['package_individual_figures']=89;metadata['package_PDF_files']=91
    (args.study/'metadata.json').write_text(json.dumps(metadata,indent=2))
    target=args.study.with_suffix('.zip')
    if target.exists():raise FileExistsError('Package already exists')
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(args.study.iterdir()):
            if path.is_file():archive.write(path,arcname=str(Path(args.study.name)/path.name))
    print(json.dumps(dict(zip=str(target),bytes=target.stat().st_size,figures=89,PDFs=91,source_sha256=digest)),flush=True)


if __name__=='__main__':main()
