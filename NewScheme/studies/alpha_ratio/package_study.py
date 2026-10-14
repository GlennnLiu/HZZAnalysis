"""Package only this study's figures, numerical results and source snapshots."""
import argparse
from pathlib import Path
import zipfile


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis', type=Path, default=Path('/cms/user/liugeliang/CMS/HZZAnalysis'))
    args = p.parse_args()
    scripts = args.analysis/'FitTemplate/AlphaRatio_examination'
    study = scripts/'plots/zh_subtraction_2018_broadCR'
    fine = scripts/'plots/zh_subtraction_2018_finecheck'
    reduced = scripts/'zhbb_reduced_2018'
    archive = study.with_suffix('.zip')
    if archive.exists(): raise FileExistsError('Refusing to overwrite archive')
    files = [(f, Path(study.name)/f.name) for f in study.iterdir() if f.is_file()]
    files += [(f, Path(study.name)/'finecheck'/f.name) for f in fine.iterdir() if f.is_file()]
    files += [(f, Path(study.name)/'zhbb_reduced'/f.name) for f in reduced.iterdir() if f.is_file()]
    files += [(scripts/name, Path(study.name)/'scripts'/name) for name in
              ['prepare_zh_reduced.py', 'compare_zh_subtraction.py', 'test_zh_subtraction.py',
               'validate_study.py', 'package_study.py']]
    files += [(args.analysis/'FitTemplate/make_reduced_trees.py', Path(study.name)/'provenance/make_reduced_trees.py'),
              (args.analysis/'cards/config_UL2018.yml', Path(study.name)/'provenance/config_UL2018.yml'),
              (Path('/cefs/higgs/zhuhao/UnbinFit/cards/median/NetSF_signal_2018.yml'), Path(study.name)/'provenance/NetSF_signal_2018.yml')]
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for source, destination in files: z.write(source, str(destination))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
    print('Saved %s: %d files, %.2f MB' % (archive, len(files), archive.stat().st_size/1e6), flush=True)


if __name__ == '__main__': main()
