# ParticleNet effect on alpha — 2018

DY only, separate ee and mumu in each of six topology/tag categories. Source:
the finalized `FitTemplate/RedTrees_20261009_pnet` trees. Both alpha versions
are freshly computed on identical selections and edges; the comparison does not
mix old trees/alphas with new ones. Original alpha ROOT files are read for their
bin edges only and are not overwritten or redeployed.

`alpha_SR = DY_SR / DY_fullCR`.

`alpha_VR = DY_VR / DY_(CR excluding VR)`.

CR retains 105–135 GeV. Selection is m2l2q >=300 GeV and 40 <= m2q <180 GeV;
SR is strict 70–105, VR strict 135–165. VR exclusion is performed event-by-event
on the CR tree, not by subtracting a potentially differently selected VR tree.
No Data, TT or VZ participates in alpha itself. No signal-region observations
are read. DY weights differ only by ParticleNet; signed weights are retained.

Use each channel/category's stored nominal SR or VR alpha edges. The first bin
is intersected with the search domain (lower bound 300); all subsequent original
edges remain unchanged. No bin merging or common ee/mumu rebinning. Plots stop
at each original histogram's finite upper edge. The integrated summary includes
all m2l2q >=300, including beyond the plotted range. Explicit total sums are
retained in the histogram archive.

Separate PNG/PDF for 24 individual comparisons; four six-category overviews
(ee/mumu, SR/VR) and two twelve-page PDF bundles. `summary.csv` contains integrated
results and maximum bin shifts in statistically populated bins; `bins.csv`
contains every bin. `histograms.npz` stores input yields/variances/covariance,
edges, alpha results and masks. `input_audit.json` records original alpha file
hashes. Code and tests are retained alongside the plots.

## Correlations and interpretation

Let A0=T0/C0 and A1=T1/C1. T is the target, C the disjoint source; subscripts
denote ParticleNet off/on. Same-event weight products provide Cov(T0,T1) and
Cov(C0,C1). First-order weighted-MC statistical propagation includes these in
Cov(A0,A1), the variance of A1/A0 and the variance of A1−A0. Do not add errors
of the two correlated alphas as if they were independent.

The change is exactly `A1/A0 = (T1/T0)/(C1/C0)`: it depends on the ratio of
weighted-average SFs, not the absolute normalization shift of DY. A common
constant SF cancels identically. SF/flavour/pT mixtures can nevertheless vary
with mass and region and generate a residual bin-dependent change.

Only positive target/source yields for both choices are plotted as ratios;
invalid bins remain in the archive. Populated-bin summaries use DY Neff >=5
in both regions and both choices, without dropping sparse bins from plots.
Only MC statistical uncertainty is included; calibration systematics and
modeling uncertainties are not. These are diagnostic nominal alpha ratios.

Archive alpha rows: A0, A1, Var(A0), Var(A1), Cov(A0,A1), A1/A0,
Var(A1/A0), Var(A1−A0). Input rows: yield0, yield1, variance0, variance1,
paired covariance. No combined-channel alpha is applied to any prediction.
The `combined_integrated_DY_ratio` summary is solely a ratio of total ee+mumu
DY yields for explaining normalization cancellation, not a replacement for
independent channel transfers.

## Reproduction

```bash
cd /cms/user/liugeliang/CMS/HZZAnalysis/Plot/alpha_pnet_comparison_2018
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python -m unittest test_alpha_pnet -v
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python check_alpha_pnet.py \
 --input /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/RedTrees_20261009_pnet \
 --alpha-dir /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/AlphaRatio \
 --output /path/to/fresh/plots
```
