"""Submit the eight authorized signal-reduction jobs once, after augmentation."""
import json
import re
import stat
from pathlib import Path
import subprocess

here=Path(__file__).resolve().parent
manifest=Path('/cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/SignalWithDeepJet_20261012/augmentation_manifest.json')
records=json.loads(manifest.read_text()); assert len(records)==24 and all(r['status']=='passed' for r in records)
target=here/'submitted_jobs.json'
if target.exists(): raise FileExistsError('Refusing duplicate submission')
jobs=[]
for era in ('2016preAPV','2016postAPV','2017','2018'):
    for cat in ('resolved','merged'):
        script=here/f'run_{era}_{cat}.sh'
        script.chmod(script.stat().st_mode|stat.S_IXUSR)
        subprocess.run(['bash','-n',str(script)],check=True)
        result=subprocess.run(['/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/hep_sub',str(script),'-g','cms','-mem','6000',
                               '-o',str(here/f'{era}_{cat}.out'),'-e',str(here/f'{era}_{cat}.err')],cwd=here,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        cluster=re.search(r'cluster\s+(\d+)',result.stdout)
        jobs.append(dict(era=era,category=cat,script=str(script),returncode=result.returncode,cluster_id=int(cluster.group(1)) if cluster else None,submission_output=result.stdout))
        target.write_text(json.dumps(jobs,indent=2))
        print(era,cat,result.stdout,flush=True)
        if result.returncode or not cluster or 'ERROR:' in result.stdout: raise RuntimeError('Submission failed; inspect submitted_jobs.json before any retry')
