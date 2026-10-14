"""Collect historical code/configuration only for the NewScheme relocation audit."""
import hashlib,json,shutil,tarfile
from pathlib import Path
BASE=Path('/cms/user/liugeliang/CMS/HZZAnalysis')
STAGE=Path('/scratchfs2/cms/liugeliang/newscheme_migration_20261014')
GROUPS={
 'production/background':'FitTemplate/jobs/reduced_20261009_pnet',
 'production/signals':'FitTemplate/jobs/signal_reduced_20261012',
 'plots/z_kinematics':'Plot/run2_z_kinematics_PNet',
 'plots/discriminants':'Plot/run2_discriminants_PNet',
 'plots/pnet':'Plot/run2_PNet_scores',
 'plots/pnet_fullrange':'Plot/run2_PNet_fullrange',
 'plots/dy_pnet':'Plot/dy_pnet_comparison_allEras',
 'plots/sr_pnet':'Plot/sr_pnet_datamc_2018',
 'plots/alpha_pnet':'Plot/alpha_pnet_comparison_2018',
 'plots/alpha_closure':'Plot/alpha_closure_pnet_2018',
 'plots/signal_overlay':'Plot/signal_overlay_Run2_20261013',
 'studies/alpha_ratio':'FitTemplate/AlphaRatio_examination',
}
def main():
    root=STAGE/'NewScheme'; root.mkdir(parents=True,exist_ok=True)
    records=[]
    for destination,original in GROUPS.items():
        source=BASE/original
        for f in sorted(source.iterdir()):
            if not f.is_file() or f.suffix not in ('.py','.md','.sh'): continue
            target=root/destination/f.name; target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(f,target)
            records.append(dict(old=str(f),new=str(target.relative_to(root)),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
        if destination.startswith('production/'):
            for asset in ('cards','sf'):
                for f in sorted((source/asset).rglob('*')):
                    if not f.is_file() or f.suffix not in ('.yml','.yaml','.json','.gz'): continue
                    target=root/destination/asset/f.relative_to(source/asset); target.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(f,target)
                    records.append(dict(old=str(f),new=str(target.relative_to(root)),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
            for name in ('submitted_jobs.json',):
                f=source/name
                if f.is_file(): shutil.copy2(f,root/destination/name)
    f=Path('/scratchfs2/cms/liugeliang/vr_no_mass_20261013/redraw_without_mass.py')
    target=root/'studies/alpha_ratio/redraw_without_mass.py'; shutil.copy2(f,target)
    records.append(dict(old=str(f),new=str(target.relative_to(root)),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    (root/'SOURCE_MANIFEST.json').write_text(json.dumps(records,indent=2))
    with tarfile.open(STAGE/'NewScheme_sources.tar.gz','w:gz') as t: t.add(root,arcname='NewScheme')
    print('Collected',len(records),'code/config files; no ROOT, caches, plots or job logs copied.')
if __name__=='__main__': main()
