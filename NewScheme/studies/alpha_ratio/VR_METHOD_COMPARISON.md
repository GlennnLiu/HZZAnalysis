# Five-method VR comparison, 2018

Compare direct nominal DY, direct pT-reweighted DY, direct mass-reweighted DY,
nominal alpha and pT-reweighted alpha to target signed Data-TTplusWW-VZ.
Independent resolved/merged x btag/untag/vbftag x ee/mumu corrections; combined
sums independently corrected/transferred channels. No ZH subtraction or
pT-plus-mass combined correction. Original mass comparison bins preserved.

Inputs: `plots/pt_reweight_2018_channel` and
`plots/mass_reweight_2018_channel_fineAlpha`. Direct VR spectra and full joint
covariance were already computed by these studies but not shown in the SR
comparison. This script validates them and recomputes all 180 VR chi2 values
from the saved residual covariance and populated-bin mask. Joint ratio errors
are independently checked from observation, prediction and residual variances.
Eight common input/nominal/covariance fields match in all 36 category/channel/
calibration combinations, and source hashes match. No fits, trees or skims
are regenerated. Missing modeling/normalization/adaptive-bin systematics remain.

Two comparison sets are kept separate:

- `inclVR`: corrections calibrated in C+VR, matching the SR primary study.
  Direct reweighted VR and pT-reweighted alpha use VR-derived factor values,
  so they are partly in-sample. Joint fit/observation covariance is included,
  but this does not make the test an independent validation.
- `heldoutVR`: fitted factor values use C only; VR observations are excluded
  from their fits. Prefer this set for a transferability cross-check. The mass
  correction boundaries were previously chosen using both calibration sets,
  including primary VR statistics; this is not a fully blind adaptive-grid
  cross-validation, and its bin-choice uncertainty is not included.

In both sets alpha uses DY(VR)/DY(C) applied to signed Z(C), not C+VR as its
source. C is accepted hadronic mass 40-180 excluding SR70-105 and VR135-165,
retaining 105-135. The entire search-mass spectrum >=300 is included. No
signal-free total-mass interval is introduced. Same-bin mass factors cancel
in alpha exactly, so mass-reweighted alpha is omitted. CR-only direct mass
reweighting equals nominal alpha on matched correction bins; finer comparison
bins can differ due to coarse-factor smoothing.

42 CMS-style figures, each PDF and PNG: 18 category/channel plus three
six-category overviews per calibration. Two 21-page bundles. Five curves with
prediction statistical bands and joint observation/prediction ratio errors.
Methods share source data and MC; bands cannot be quadrature-added to test
their differences. Statistics use full residual covariance, retain absolute
normalization, and apply the existing populated-bin masks (Data>=5, nominal
DY Neff>=5). Denominator in chi2/rank is residual covariance rank, not a count
of reweighting fit parameters subtracted from VR bins. Tail bins remain visible
in spectra; no high-mass hatching or exclusion. Last bin includes overflow.

`summary.csv` contains yields, chi2/rank and p-values in every channel/category.
`<channel>_VR_chi2_<inclVR|heldoutVR>.md` gives compact five-method tables.
`results.npz` stores all plotted values, full prediction/residual covariance,
observations and joint ratio errors. Metadata records input archives/hashes.
No final-likelihood or deployment claim.

TTplusWW=TTTo2L2Nu+WWTo2L2Nu+tZq; VZ=WZTo2Q2L+ZZTo2Q2L, preserving the
producer's existing Data-TT-VZ shorthand.

Reproduce from AlphaRatio_examination with fresh output:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python plot_vr_methods.py \
  --pt plots/pt_reweight_2018_channel \
  --mass plots/mass_reweight_2018_channel_fineAlpha \
  --output plots/VR_method_comparison_2018_NEW
```

Original result folders are required for reproduction; their full caches are
not duplicated into this plotting package. Current files and ROOT producers
are unchanged by this presentation/validation study.
