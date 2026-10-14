# Run-2 ParticleNet scores: DY SF on/off

AK8-candidate scores particleNetZvsQCD and particleNetZbbvslight, using the merged selection. Resolved dijets have no corresponding ParticleNet score; a stored AK8 score in a resolved event would describe a different object and is not plotted as a dijet tagger.

Combine all four eras, all tags and ee+mumu. Retain original selections and normalizations, m2l2q>=300 and 40<=mass2q<180. Label luminosity as 138 fb^-1. CR includes 105–135 and VR, labelled simply CR; do not add VR trees separately. SR Data is never opened. Direct MC is used, not alpha estimation.

The VZ/TT/DY colored stack uses nominal weights including DY ParticleNet SF. Green outline shows its total; red dashed outline shows the total with only DY switched to weight_no_pnet. Other MC stays identical in both predictions. In CR, black ratio markers use the SF-on denominator; red open squares use the SF-off denominator. Gray bands show nominal MC statistics (sumw2), not systematics. Red dashed ratio-band boundaries show SF-off MC statistics. Both predictions share the same MC events and observed data and are correlated; these bands are not independent uncertainty estimates for their difference.

20 bins from 0–1, with exact endpoints included. Nonfinite/out-of-range scores are excluded and audited. These are selected merged-event distributions, not preselection tagger-efficiency plots: the existing ZvsQCD working-point requirement is retained.

Four separate PNG/PDF comparisons and two overview PNG/PDF pairs, with a four-page PDF bundle. Script, histogram cache, input provenance and per-input audit are retained. Validation checks full weighted yields/variances, no SR Data, and eventwise weight=weight_no_pnet*pnet_sf for DY.

```bash
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python plot_run2_pnet.py --input /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/RedTrees_20261009_pnet --output .
```

For plot-only updates add `--histograms histograms.npz`. Previous plotting packages are untouched.
