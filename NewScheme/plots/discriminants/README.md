# Run-2 discriminants

All four eras, ee+mumu and all exclusive tags combined, separately resolved/merged. Full nominal MC weights include DY ParticleNet SFs. Direct DY MC, not alpha estimation. Original era normalization is retained; luminosity label is 138 fb^-1. Selections match the Run-2 Z kinematic package (m2l2q>=300, 40<=mass2q<180). CR includes 105–135 and VR without adding VR trees; the label is simply CR. SR Data trees are never opened.

ANATree/src/TreeLoop.cc computes jjVBF for >=2 additional jets, and jVBF only for exactly 1 additional jet (else-if branch). Reduced-tree aliases select KD_jjVBF/KD_jVBF in resolved and KD_JJVBF/KD_JVBF in merged; kd selects KD_Zjj/KD_ZJ. Only finite values in [0,1] are plotted, including exact endpoints. Negative unset sentinels and failed MELA values are excluded, not folded into zero. The reduced trees do not retain the associated-jet count, so the validity mask uses stored discriminant values rather than inventing a jet-count cut. The audit reports exclusion counts, valid yields, variances, and overlap of the two VBF validity masks for each input.

20 equal bins, 0–1. MC bands are sqrt(sumw2), statistical only. CR Data errors are Garwood Poisson intervals, with Data/MC panels. Twelve individual PNG/PDF pairs, four three-variable overview pairs, and one 12-page PDF bundle. Script, histograms, input metadata and validity_audit.csv are retained for reproduction and fast cached rerendering.

```bash
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python plot_run2_kd.py --input /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/RedTrees_20261009_pnet --output .
```

For layout-only updates add `--histograms histograms.npz`. Previous kinematic plots are unchanged.
