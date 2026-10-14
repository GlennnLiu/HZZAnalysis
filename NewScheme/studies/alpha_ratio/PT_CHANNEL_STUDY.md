# 2018 category- and lepton-channel-dependent DY pT reweighting

Created 2026-10-07. Output: `plots/pt_reweight_2018_channel/`, mirrored locally
and on LXLOGIN. This is a separate diagnostic, not a production correction.
All previous outputs remain unchanged. No new source scan or ROOT production.

## Definition

For each resolved/merged x btag/untag/vbftag x ee/mumu, separately derive
`r_i=(Data-TTplusWW-VZ)_i/DY_i` in hadronic-pT bins. No pooling of tags and
no combined-channel factor. Each combined-channel prediction is the sum of
its independently corrected and transferred ee and mumu predictions.

The existing producer selection, nominal weights and sample grouping are
preserved. Actual variables are pt2jet (resolved) and ptmerged (merged).
ZH is not subtracted. TTplusWW includes TTTo2L2Nu, WWTo2L2Nu and tZq; VZ
includes WZTo2Q2L and ZZTo2Q2L. Signed yields and variances are preserved.
The input cache is copied into the package with its source audit.

Hadronic SR is 70<q<105 and VR is 135<q<165 within 40<q<180. C is outside
both. Primary calibration uses C+VR; the cross-check (`heldoutVR`) uses C.
The primary SR alpha source is C+VR. The VR alpha source is C and is always
disjoint from VR; primary VR nevertheless enters calibration, so it is partly
in-sample. The `VR_alpha_heldoutVR` plots exclude VR from both calibration
and alpha source. SR alpha using C alone is also saved numerically.
All m(2l2q)>=300 remains in the search domain; there is no total-mass CR.

Adjacent pT bins are merged within each category and channel until signed Z
and DY are positive with Neff>=25 each. No target-SR data enters fitting or
pT-bin choice. ee and mumu can have different pT grids. Combined pT plots
use their common coarsening, without changing either channel's factors.
The earlier study's mass grid is retained, with additional merging only if
required by nominal/corrected DY positivity per channel and region, never
observed target yields. Nominal integrated yields/errors on unchanged grids
are checked against the earlier combined-channel-calibration study. Any
additional merging is recorded explicitly in metadata; nominal alpha yields
can depend on the mass grid and must not be attributed to the correction.

## Statistics and validation

The existing analytic first-order covariance propagates signed-data,
subtraction-MC and DY statistics through factors, weighted DY and alpha.
Calibration/data reuse is retained, including the partly in-sample VR test.
Independent channel corrections have disjoint primitives; combined predictions
sum channel Estimates before covariance propagation. Calibration before/after
closure is exact in each lepton channel and category, not just ee+mumu.
The fitted sample's ratio error vanishes with the shared-data covariance;
this is a construction check, not independent transfer validation.

Tests cover per-channel closure and isolation, separate-channel covariance,
no target-SR/held-out-VR data in prediction, finite-difference gradients,
common-grid coarsening and constant-factor cancellation. Runtime checks assert
all 36 calibration spectra close in values and Jacobians and all 36 SR
pT/mass projections have identical integrated values and derivatives.

Statistical-only Gaussian chi2 uses the full residual covariance and populated
bins (Data count>=5 and DY Neff>=5). All tail bins remain in top distributions
and tables. These are not final likelihood tests. Missing uncertainties:
lumi, tags, subtraction normalization/shape, DY modeling, normalization-induced
MC correlations and adaptive-bin-choice uncertainty. No new physics inference
or production deployment is implied by this replotting request.

## Numerical outcome

All 13 regression tests passed. All six categories retained the earlier mass
grid, so nominal mass yields/errors are unchanged. All 36 calibration spectra
and 36 SR pT/mass projection checks passed.

For untagged resolved SR, direct yield closure (Data-TTplusWW-VZ)/DY changes
from 0.89784 to 0.99652 in ee and from 0.97619 to 0.99416 in mumu. The
combined direct closure is 0.99517. Combined SR alpha closure is 0.99819 and
VR alpha closure is 0.91917. Their mass-distribution chi2/rank is 119.61/10,
62.54/10 and 159.74/10 respectively, statistical only.

For untagged merged SR, direct closure changes from 1.17637 to 0.98146 (ee)
and 1.23248 to 0.97564 (mumu). Combined direct/SR-alpha/VR-alpha closure is
0.97811/0.97844/1.21209, with chi2/rank 1.29/4, 1.21/4 and 6.47/4.
The combined untagged conclusion barely changes from the previous shared-
lepton-channel, category-specific fit. This is a channel-dependent remake,
not a demonstrated solution to alpha nonclosure. In primary VR results the
target is partly in-sample; held-out-VR plots are provided separately.

## Plots

147 individual CMS-style figures, each in PNG and PDF; no high-mass hatching.
Orange: nominal DY or nominal alpha. Dashed red: channel/category-specific
reweighted DY or alpha. Gray: normalization-only, independently per channel,
in direct SR plots. Black: Data-TTplusWW-VZ in the compared region.
There is no shared-tag/category-inclusive correction curve.

- `combined_SR_direct_summary.pdf`, `combined_SR_alpha_summary.pdf`,
  `combined_VR_alpha_summary.pdf`: six-category overviews, summing channels.
- `ee_..._summary.pdf` and `mumu_..._summary.pdf`: corresponding channel overviews.
- `<topology>_<tag>_<channel>_SR_direct.pdf`, `_SR_alpha.pdf`, `_VR_alpha.pdf`:
  individual mass distributions, with channel ee/mumu/combined.
- `<topology>_<tag>_<channel>_VR_alpha_heldoutVR.pdf`: independent VR check.
- `<topology>_<tag>_<channel>_SR_direct_pt.pdf`: before/after SR pT distributions.
- `<topology>_<tag>_<ee|mumu>_pt_scale_factors.pdf`: fitted SFs, C+VR versus C.
- `<topology>_<tag>_<channel>_<inclVR|heldoutVR>_direct_pt.pdf`: calibration pT
  distributions; category- and channel-specific exact closure by construction.
- `all_pt_reweight_plots.pdf`: 93-page mass-comparison and SF bundle.
- `all_SR_pt_reweight.pdf`: 18-page SR pT bundle.
- `all_CR_pt_reweight.pdf`: 36-page calibration pT bundle.
- `summary.csv`, `bins.csv`, `pt_factors.csv`, `results.npz`: numerical results.
- `metadata.json`: definitions, checks, script/input hashes, adaptive binning.

## Reproduction

From AlphaRatio_examination on LXLOGIN, with a fresh output path:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PYTHON=/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python
$PYTHON -m unittest test_pt_reweight test_pt_reweight_channel -v
$PYTHON pt_reweight_channel.py --source plots/pt_reweight_2018_fullCR \
  --output plots/pt_reweight_2018_channel_NEW
$PYTHON plot_pt_reweight_channel.py --input plots/pt_reweight_2018_channel_NEW
$PYTHON finalize_pt_reweight_channel.py --study plots/pt_reweight_2018_channel_NEW
```

The plotting step copies code/README and creates the sibling ZIP. Finalization
verifies completeness, reruns tests into `tests.log`, saves `validation.json`,
and refreshes the ZIP. The package
includes `baseline_summary.csv` and `baseline_mass_edges.npz` from the earlier
study and the unchanged source cache. To reproduce using only this package,
run its scripts with `--source .`
and a fresh external output path, then plot that new path. No original ROOT
inputs are needed for the cached analysis.
