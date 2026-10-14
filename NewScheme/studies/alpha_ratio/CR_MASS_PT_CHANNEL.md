# Calibration-region mass shapes after hadronic-pT reweighting

2018, channel- and category-specific corrections from `pt_reweight_2018_channel`.
Output `plots/pt_reweight_2018_channel_CR_mass/`, local and LXLOGIN mirrors.
No changes to fitted factors, existing mass grids, ROOT inputs or production.

The primary `inclVR` figures compare nominal DY, pT-reweighted DY and a
normalization-only baseline to Data-TTplusWW-VZ in the calibration region C+VR.
The `heldoutVR` figures repeat the comparison in C using the C-only factors.
No ZH subtraction. Actual mass is m2l2q resolved or m2lJ merged. Entire spectrum
above 300 is retained, including the overflow in the last labelled bin.

Both calibrations retain hadronic masses 105-135. C excludes hadronic SR
70-105 and VR 135-165 from accepted 40-180; C+VR excludes only SR.
Correction factors are independently derived in each topology/tag/ee or mumu.
Combined-channel curves sum independent corrections; there is no shared-tag
or combined-channel calibration. Individual and overview plots cover all six
categories and ee/mumu/combined; every figure has a separate PNG and PDF.

Orange is nominal DY, dashed red pT-reweighted DY, dashed gray normalization
only, black Data-TTplusWW-VZ. Lower panels compare observed/predicted with
joint statistical errors, including the dependence on the same calibration
data. Gray triangles provide the normalization-only baseline in the lower panel.

Exact pT closure does NOT enforce mass-shape closure. Integrated yields of
both pT-corrected and normalization-only predictions must agree with calibration
Z by construction; all 36 value AND Jacobian assertions are checked. A remaining
mass shape residual probes conditional correlations within pT bins, but is still
an in-sample diagnostic, not independent proof of SR transferability.

`summary.csv` reports yields and full-covariance chi2/rank. Fits introduce
constraints in the mass residual covariance: use its eigen-rank, not simply
the number of bins, and do not quadrature-add independent data/fitted-MC
uncertainties. Normalization-only already removes the total-yield mismatch;
compare its mass shape to pT reweighting to identify any additional shape gain.
Tests use populated bins (Data count>=5, nominal DY Neff>=5); all tail bins
remain in the plots/tables. Statistical-only Gaussian checks, not a final likelihood;
subtraction/tag/DY-modeling systematics and MC normalization covariance omitted.

## Outcome for untagged C+VR calibration

Combined ee+mumu mass chi2/rank:

| Topology | Nominal DY | Normalization only | pT reweighted |
|---|---:|---:|---:|
| Resolved | 384.62/10 | 37.58/9 | 70.78/9 |
| Merged | 79.84/4 | 1.45/3 | 0.62/3 |

The total observed/predicted yield is 1 for both corrected options by
construction (nominal 0.94415 resolved and 1.23893 merged). Resolved mass
shape remains inconsistent on a statistical-only test; the pT correction is
actually worse than the normalization-only baseline. This is also true
within each channel: ee norm-only 36.42/9 versus pT 57.51/9; mumu norm-only
11.48/9 versus pT 23.97/9. Merged is compatible at this coarse binning and
statistical precision after either correction. This does not establish that
the correction transfers to SR or uniquely identify the underlying cause.

## Files

- `all_CR_mass_pt_reweight.pdf`: 42-page bundle (21 primary plus 21 C-only).
- `combined_inclVR_mass_summary.pdf`: all six categories, primary C+VR, ee+mumu.
- `ee_inclVR_mass_summary.pdf`, `mumu_inclVR_mass_summary.pdf`: channel overviews.
- `<topology>_<tag>_<channel>_<inclVR|heldoutVR>_mass.pdf`: individual distributions.
- `summary.csv`, `bins.csv`, `results.npz`: numerical outputs and covariances.
- `validation.json`: unchanged-input/factor and integrated-closure checks.
- Source cache, reference factors/metadata, helpers and script are included.

Reproduction, from AlphaRatio_examination, using a fresh output folder:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python \
  plot_cr_mass_pt_reweight_channel.py --study plots/pt_reweight_2018_channel \
  --output plots/pt_reweight_2018_channel_CR_mass_NEW
```

From within this package, use `--study .` with an external fresh output path.
