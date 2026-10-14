# Signal production handoff

Requested ggF and VBF at 500/1000/3000 GeV for 2016preAPV, 2016postAPV, 2017 and 2018. All 24 mass-specific Zhuhao skims exist, including usable VBF3000 alternatives for the two empty Guojl skims. Original ntuples/YAMLs/shared helpers/background reduced trees remain untouched.

2026-10-12: all-era nominal/up/down two-tag DeepJet calibration regression passed against corrected combined signal files. Signal-specific efficiencies and event-order-preserving pass/fail weights are calculated separately for each mass point. Efficiencies/counts and calibration provenance are saved. Full cloned signal skims preserve entry order, unrelated branches and normalization, adding w_deepjet, w_deepjet_up, w_deepjet_dn, w_deepjet_valid. Domain-invalid weights are explicit -999 sentinels, never multiplied into selected VBF/merged weights.

Completed and validated **24 augmented files**, **3547095 source entries**, **355668 DeepJet-valid entries**. Nominal weight extrema 0.4208441719–2.4724916857; no negative nominal or calibration-variation weights. Output: `/cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/SignalWithDeepJet_20261012/`. Manifest and all per-sample reports are saved there. The original file sizes/mtimes, source IDs and key physics arrays were checked against the cloned outputs.

The ggF500 2018 pilot passed resolved and merged reduction: 18 files per topology, fully read and validated. Full reduction uses frozen background-production cut aliases, retaining 105–135 in full CR and VR=135–165. Merged signal gets era/flavour/pT/WP ParticleNet factors, resolved btag/untag gets nominal DeepJet. Save weight_base, weight_no_pnet, weight_no_deepjet, full weight, both PNet factors, raw DeepJet variations/validity, nominal/generated mass, SumWeight and source entry. Signal card cross sections are 1 pb benchmark normalization.

Eight full-production jobs are submitted. The first submission attempt was rejected for missing script executable permission and queued no jobs; its audit is retained as `submitted_jobs.failed_permissions.json`. Executable permissions were added to the new scripts only, and the successful submission script now checks for actual cluster IDs (the hep_sub wrapper can return exit code zero on errors).

| Era | Resolved job | Merged job |
|---|---|---|
| 2016preAPV | 82173109 | 82173110 |
| 2016postAPV | 82173111 | 82173112 |
| 2017 | 82173113 | 82173114 |
| 2018 | 82173115 | 82173116 |

Latest status: **all eight queued/idle**, no full-production validation reports yet. Expected **864 reduced ROOT files**, 108 per job, under `/cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/RedTrees_signal_20261012/`. Each producer independently rereads and validates every file before marking that era/topology done. Do not treat submission, or the pilot, as full production completion.

When all eight `.done` markers exist, run:

```bash
source /scratchfs2/cms/liugeliang/miniconda3/bin/activate unbin
cd /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/jobs/signal_reduced_20261012
python validate_complete.py --augmented /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/SignalWithDeepJet_20261012 --reduced /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/RedTrees_signal_20261012
```

This saves `complete_validation.json` only after all 24 augmented and 864 reduced products pass, including original-input and frozen-YAML preservation checks. If a batch job fails, inspect its `<era>_<topology>.out/.err` before retrying; producers refuse overwriting existing final or partial ROOT files. Do not resubmit the successful jobs wholesale.
