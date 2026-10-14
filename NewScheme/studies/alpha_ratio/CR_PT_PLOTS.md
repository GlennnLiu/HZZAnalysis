# 2018 calibration-region pT before/after plots

Saved 2026-10-07. This supplements `pt_reweight_2018_fullCR`; no correction,
ROOT input, selection, weight or production alpha is changed.

`all_CR_pt_reweight.pdf` contains 48 pages. Every page also has a separate
CMS-style PDF and PNG. The first 24 use the main C+VR calibration (`inclVR`);
the last 24 use the C-only cross-check (`CRonly`). Each set contains resolved
and merged, three separate tags and their pooled sum, combined ee+mumu and
each individual channel. Actual variables are pt2jet and ptmerged.

Nominal DY is orange. Blue is the correction fitted with all three tags pooled.
For separate-tag pages, dashed red uses that tag's own correction. Black points
are Data-TTplusWW-VZ (the existing producer groups); ZH is not subtracted.
Each page uses that fitted sample's adaptive pT bins, not a finer arbitrary grid.
Sparse adjacent bins were merged using the original Neff>=25 rule.

The correction is r_i = Z_i / DY_i in the selected calibration sample, with
Z=Data-TTplusWW-VZ and ee+mumu combined. Thus corrected DY_i = Z_i exactly
for that combined fitted sample. This is enforced by construction, not a test
of the correction's transferability to SR. A common correction need not close
for every tag or individual channel. Likewise C+VR closure does not enforce
C-only or VR-only closure separately.

Statistical bands and ratio errors use the existing analytic first-order
Jacobian, retaining calibration/observation covariance. For the fitted sample
both yields and Jacobians agree; its corrected observed/predicted ratio is
one with vanishing propagated statistical error. Treating fitted DY and the
same calibration data as independent would double count uncertainty.
Systematic uncertainties are not included. Open-ended last bins are labelled
with their lower bound; there is no high-pT hatching.

`cr_pt_summary.csv` contains yields/closure; `cr_pt_bins.csv` contains bin values.
`cr_pt_distributions.npz` stores sumw/variance rows and covariance matrices.
`validation.json` records the input hash and all 16 combined-sample assertions.

Reproduce from the existing source cache in AlphaRatio_examination:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python plot_cr_pt_reweight.py \
  --source plots/pt_reweight_2018_fullCR \
  --output plots/pt_reweight_2018_CR_closure_NEW
```

The output path must be fresh. The script copies its helpers and this README
into the output folder and creates a sibling ZIP. The reusable source cache
remains in `pt_reweight_2018_fullCR/source_histograms.npz`.
