# 2018 CR mass reweighting: direct and alpha closure

2026-10-07. Primary output `plots/mass_reweight_2018_channel_fineAlpha/`, local
and LXLOGIN copies. `mass_reweight_2018_channel` is the earlier coarse-alpha-
grid diagnostic; use the fineAlpha version for comparison with previous plots.
Independent resolved/merged x btag/untag/vbftag x ee/mumu corrections. No
category-inclusive fit, no pT reweighting, no ZH subtraction, no production
or ROOT changes. The original nominal-DY source cache/producer weights and
selection are used, not the already pT-reweighted distributions.

## Construction and identities

In each mass bin derive r_i = Z_cal,i / DY_cal,i, with Z=Data-TTplusWW-VZ.
Apply the same r_i to nominal DY in SR and VR. Primary calibration is C+VR;
the held-out-VR check fits C alone. SR is hadronic 70-105, VR 135-165, accepted
40-180; C excludes both. Both retain 105-135. Primary SR alpha uses C+VR;
VR alpha uses C, disjoint from VR. Primary VR enters calibration, so a direct
VR comparison is partly in-sample. In the alpha prediction the mass-factor
dependence cancels completely: calibration adds no VR-data dependence there.
The whole search mass spectrum >=300 is retained.

The reweighting variable is m2l2q (resolved) or m2lJ (merged). Its correction grid begins
with the previous mass edges and is merged until Z and DY have Neff>=25 in
each calibration/channel and target DY is positive. The most sparse channel
sets shared bin boundaries, not shared correction values. Each channel/category
is fitted independently. No SR observations choose bins or correction values.
VR observations enter only the primary calibration, not the C-only cross-check.
The original alpha/comparison grid is preserved without any rebinning. Each
coarse mass correction is expanded as a constant over its constituent alpha
bins; no alpha bin crosses a correction boundary. Nominal predictions and
their integrated errors are checked against the previous channel-dependent
study. Thus before/after alpha comparison is not confounded by coarser bins.
Adaptive-bin uncertainty and signal leakage into sidebands are not established.

Therefore, in each channel and bin:

```
alpha_rw = (r * DY_target) / (r * DY_source) = alpha_nominal
```

If the correction is derived from that source, on its correction grid:

```
direct_DY_rw = (Z_source/DY_source) * DY_target
             = alpha_nominal * Z_source
```

Both equalities hold for values AND propagated derivatives/covariances, not
just central yields. They are verified at runtime. Combined plots sum separately
calibrated/transferred channels, not a combined fit or combined alpha.

Consequently, direct mass-reweighted SR DY is the usual nominal alpha estimate
when both use the correction grid. Where factors are smoothed into coarser bins,
the displayed fine-bin direct prediction can differ from fine-bin alpha: that
is an explicit smoothing assumption. Recomputing alpha from the mass-reweighted
DY still cannot improve closure on the preserved alpha grid. Any apparent improvement
in alpha from mismatched/coarse rebinning must be audited as a bin-averaging
effect, not claimed as a new transferable mass correction. Smooth corrections
inside finite bins need the analogous averaging audit.

## Statistics and validation

Signed data/background subtraction and nominal DY sumw/sumw2 primitives feed
the existing joint first-order Jacobian covariance. Fit/data reuse is retained.
Fitted CR mass distributions close in their correction bins with vanishing
joint ratio error; the fine-bin plots need not close within coarse correction
bins. Exact fitted-grid closure is imposed, not independent validation. Full residual covariance is used
for statistical-only Gaussian chi2 on populated bins (Data count>=5 and DY
Neff>=5). All tail bins remain plotted and tabulated. No final likelihood claim.
Lumi/tag/subtraction/DY-modeling systematics and normalization-induced MC
correlations are not included; adaptive grid choice is treated as fixed.

Before and after curves use the SAME ORIGINAL alpha/comparison grid. Sparse
categories have coarser correction bins; `metadata.json` records both sets of
edges. Do not credit correction smoothing or any coarse-alpha-grid diagnostic
changes in nominal alpha to an improved transfer factor.

15 regression tests (seven original plus eight mass-specific) check gradients,
channel isolation, no SR-data bin choice, source/target disjointness, CR closure,
nonconstant-factor alpha cancellation and the direct-equals-alpha identity.
Runtime checks cover 36 CR closures, 72 SR/VR alpha cancellations and 36 matched-
source direct/alpha identities on the correction grid, with ee/mumu/combined
in both calibrations. Fine-bin calibration yields/derivatives also close when
integrated, and nominal SR/VR predictions agree with the earlier study.
The package finalizer reruns tests and verifies all alpha covariance matrices.

## Measured outcome (untagged, combined ee+mumu)

| Test | Nominal observed/predicted | Mass-reweighted | Nominal chi2/rank | Reweighted chi2/rank |
|---|---:|---:|---:|---:|
| Resolved direct SR | 0.94088 | 0.99951 | 253.38/10 | 65.43/10 |
| Merged direct SR | 1.20793 | 0.97482 | 29.80/4 | 1.33/4 |
| Resolved alpha SR | 0.99951 | 0.99951 | 65.42/10 | 65.42/10 |
| Merged alpha SR | 0.97498 | 0.97498 | 1.20/4 | 1.20/4 |
| Resolved alpha VR | 0.91705 | 0.91705 | 167.24/10 | 167.24/10 |
| Merged alpha VR | 1.22824 | 1.22824 | 6.59/4 | 6.59/4 |

Results are statistical only, with joint observation/prediction covariance.
The corrected direct SR ratios have errors 0.00633 resolved and 0.03913 merged.
Normalization-only direct SR gives chi2/rank 127.06/10 resolved and 1.56/4 merged.
The mass correction improves absolute DY but supplies no new improvement to
the alpha estimate: its common numerator/denominator factor cancels.
Resolved direct shape still fails a stat-only compatibility check; unit integrated
yield does not imply mass-shape closure. Merged is compatible at this precision.

Correction bins, resolved/merged: btag 4/2, untag 9/3, vbftag 7/1, shared edges
but independent channel values. Merged VBF correction is normalization-only
because its calibration statistics support only one bin. All original plotted
mass grids and nominal alpha predictions remain unchanged.

## Plots and files

135 individual CMS-style figures, each PDF and PNG; no high-mass hatching.
Orange: nominal DY or nominal alpha. Red dashed: mass-reweighted DY or alpha.
Gray: normalization-only for direct comparisons. Black: Data-TTplusWW-VZ in
the compared region. Nominal and corrected alpha overlap exactly.

- `combined_SR_direct_summary.pdf`, `combined_SR_alpha_summary.pdf`,
  `combined_VR_alpha_summary.pdf`: six-category overviews (channel sum).
- `ee_..._summary.pdf`, `mumu_..._summary.pdf`: individual-channel overviews.
- `<topology>_<tag>_<channel>_SR_direct.pdf`, `_SR_alpha.pdf`, `_VR_alpha.pdf`.
- `<topology>_<tag>_<channel>_VR_alpha_heldoutVR.pdf`: C-only calibration check.
- `<topology>_<tag>_<ee|mumu>_mass_scale_factors.pdf`: fitted mass factors.
- `<topology>_<tag>_<channel>_<inclVR|heldoutVR>_CR_mass.pdf`: exact CR closure.
- `all_mass_reweight_plots.pdf`: 93-page SR/VR comparisons and factor bundle.
- `all_CR_mass_reweight.pdf`: 42-page calibration bundle, including six overviews.
- `summary.csv`, `bins.csv`, `mass_factors.csv`, `results.npz`: spectra/statistics.
- `metadata.json`, `validation.json`, `tests.log`: definitions, grids and checks.
- Original source cache/audit, baseline mass edges, scripts/helpers are included.

TTplusWW = TTTo2L2Nu+WWTo2L2Nu+tZq; VZ = WZTo2Q2L+ZZTo2Q2L, preserving the
analysis's existing grouping used in the user's Data-TT-VZ shorthand.

## Reproduction

From AlphaRatio_examination, with fresh output:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PYTHON=/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python
$PYTHON -m unittest test_pt_reweight test_mass_reweight -v
$PYTHON mass_reweight_dy.py --source plots/pt_reweight_2018_channel \
  --output plots/mass_reweight_2018_channel_fineAlpha_NEW
$PYTHON plot_mass_reweight.py --input plots/mass_reweight_2018_channel_fineAlpha_NEW
```

From within this self-contained package use `--source .` and a fresh external
output folder. Finalize without regenerating plots using `plot_mass_reweight.py
--input OUTPUT --finalize-only`. Its sibling ZIP contains the complete package.
