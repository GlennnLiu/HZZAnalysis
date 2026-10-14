"""Install reversible lxlogin aliases after migration validation, preserving originals."""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent
def main():
    p=argparse.ArgumentParser(); p.add_argument('--apply',action='store_true'); args=p.parse_args()
    assert json.loads((ROOT/'MIGRATION_VALIDATION.json').read_text())['status']=='passed'
    plan=[]
    for record in json.loads((ROOT/'SOURCE_MANIFEST.json').read_text()):
        target=ROOT/record['new']; old=Path(record['old'])
        if target.suffix!='.py' or record['new'].startswith('production/'): continue
        assert target.is_file() and ROOT in target.resolve().parents
        assert BASE in old.parent.resolve().parents or str(old.parent)=='/scratchfs2/cms/liugeliang/vr_no_mass_20261013'
        if old.is_symlink() and old.resolve()==target.resolve(): continue
        assert old.is_file() and not old.is_symlink()
        assert hashlib.sha256(old.read_bytes()).hexdigest()==record['sha256'],('Changed source',old)
        backup=old.with_name(old.name+'.pre_NewScheme_20261014'); assert not backup.exists()
        plan.append((old,target,backup))
    if not args.apply:
        print('Validated compatibility plan for',len(plan),'scripts; no originals changed.'); return
    applied=[]
    try:
        for old,target,backup in plan:
            old.rename(backup)
            try: old.symlink_to(target)
            except Exception: backup.rename(old); raise
            applied.append((old,target,backup))
    except Exception:
        for old,target,backup in reversed(applied): old.unlink(); backup.rename(old)
        raise
    report=[dict(old=str(a),canonical=str(b),preserved_original=str(c)) for a,b,c in applied]
    (ROOT/'COMPATIBILITY.json').write_text(json.dumps(report,indent=2))
    print('Installed',len(report),'compatibility aliases; originals preserved, frozen production unchanged.')
if __name__=='__main__': main()
