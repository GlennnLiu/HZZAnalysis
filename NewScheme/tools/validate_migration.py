"""Non-producing checks for relocated code, helpers and frozen assets."""
import ast,hashlib,json,os,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    env=os.environ.copy()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): env[key]='1'
    env['MPLCONFIGDIR']=str(BASE/'Plot/signal_overlay_Run2_20261013/mplconfig')
    files=list(ROOT.rglob('*.py'))
    for f in files: compile(f.read_bytes(),str(f),'exec')
    source=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text()); assets=0
    for record in source:
        target=ROOT/record['new']
        if target.suffix in ('.yml','.yaml','.json','.gz'):
            assert sha(target)==record['sha256']; assets+=1
    old=BASE/'FitTemplate/jobs/reduced_20261009_pnet/make_reduced_trees.py'
    def functions(p): return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(p.read_text()).body if isinstance(n,ast.FunctionDef)}
    before=functions(old); after=functions(ROOT/'production/background/make_reduced_trees.py')
    names=('needs_particle_net','selection_branches','input_path','particle_net_sf','mc_weights','reduced_record','apply_alpha','make_zjet')
    for name in names: assert before[name]==after[name],name
    commands=['production/background/produce_reduced_stream.py','production/background/make_reduced_trees.py',
              'production/signals/prepare_signal.py','production/signals/produce_signal_reduced.py',
              'plots/z_kinematics/plot_run2_z.py','plots/discriminants/plot_run2_kd.py','plots/pnet/plot_run2_pnet.py',
              'plots/pnet_fullrange/plot_fullrange_pnet.py','plots/dy_pnet/plot_dy_pnet.py','plots/sr_pnet/check_sr_pnet.py',
              'plots/alpha_pnet/check_alpha_pnet.py','plots/alpha_closure/check_alpha_closure.py',
              'plots/signal_overlay/prepare_histograms.py','plots/signal_overlay/render_plots.py','plots/signal_overlay/revise_layouts.py',
              'studies/alpha_ratio/pt_reweight_channel.py','studies/alpha_ratio/plot_sr_methods.py','studies/alpha_ratio/plot_vr_methods.py',
              'studies/alpha_ratio/plot_sr_zh_statband.py','studies/alpha_ratio/redraw_without_mass.py']
    results=[]
    with tempfile.TemporaryDirectory(prefix='hzz_newscheme_check_') as cwd:
        for name in commands:
            result=subprocess.run([sys.executable,str(ROOT/name),'--help'],cwd=cwd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=90)
            assert result.returncode==0,(name,result.stdout)
            results.append(name); print('CLI passed:',name,flush=True)
        for name in ('production/background','studies/alpha_ratio'):
            patterns=('test_reduced_pnet.py',) if 'background' in name else ('test_pt_reweight.py','test_pt_reweight_channel.py')
            for pattern in patterns:
                result=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(ROOT/name),'-p',pattern],cwd=cwd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=90)
                assert result.returncode==0,(pattern,result.stdout); print(result.stdout,flush=True)
    for f in ROOT.glob('production/*/*.sh'):
        subprocess.run(['bash','-n',str(f)],check=True)
    report=dict(status='passed',python_files=len(files),cli_help_checked=results,frozen_assets_identical=assets,
                unchanged_physics_functions=list(names),no_production=True,no_refit=True,no_submission=True)
    (ROOT/'MIGRATION_VALIDATION.json').write_text(json.dumps(report,indent=2))
    print('PASS: relocation syntax, imports, CLI, pure weight/covariance tests and frozen calibration hashes.',flush=True)
if __name__=='__main__': main()
