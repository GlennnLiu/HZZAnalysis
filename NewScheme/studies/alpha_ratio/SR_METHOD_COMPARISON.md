# Five-method SR comparison, 2018

Compare nominal direct DY, direct pT-reweighted DY, direct mass-reweighted DY,
nominal alpha estimation and pT-reweighted alpha estimation with signed
Data-TTplusWW-VZ in the SR. No ZH subtraction. Each topology/tag/lepton
channel is corrected independently; combined sums separately fitted ee/mumu.

Uses `plots/pt_reweight_2018_channel` and
`plots/mass_reweight_2018_channel_fineAlpha`, not the earlier shared-channel
pT study or coarse-alpha mass diagnostic. Fits and trees are not regenerated.
The script checks common primitive-source hash, original mass edges, signed
observation, observation covariance, comparison mask, nominal DY and nominal
alpha values/full covariances in all 18 topology/tag/channel combinations.

Primary corrections and nominal SR alpha source use hadronic C+VR. SR
hadronic mass is 70-105 GeV, VR 135-165, accepted 40-180; C excludes both,
retaining 105-135. Entire search mass >=300 is included. No SR observation
was used to derive corrections. Mass correction bins may be coarser than
the preserved display grid, with constant factors expanded over original bins.

Each of 18 individual figures is saved in PNG and PDF, along with three
six-category ee/mumu/combined overviews and a 21-page bundle. CMS Work in
progress style. All five curves are retained even where they nearly coincide.
Bands show prediction statistical uncertainties. Colored lower-panel points
show observation/prediction with joint statistical covariance, offset slightly
within bins for visibility. Only populated stable ratio bins are shown, matching
the earlier studies; all spectrum bins remain in the upper plot and archive.
The open-ended last bin includes overflow; no high-mass hatching or exclusion.

Integrated observed/predicted ratios are included in individual legends.
`summary.csv` preserves original closure errors, chi2/rank and p-values.
`results.npz` contains every plotted spectrum, covariance and ratio. Methods
share data/MC and are strongly correlated: their bands are not independent
uncertainties and must not be quadrature-added to test method differences.
No modeling, subtraction, lumi/tag, adaptive-bin or signal-leakage systematics
are established by this plotting comparison.

Mass-reweighted alpha is deliberately omitted: the common same-bin mass factor
cancels exactly in target/source DY, including propagated covariance. On the
correction grid direct mass reweighting equals nominal alpha transfer; on the
finer display grid explicit smoothing can cause small differences. Five curves
therefore do not represent five independent constraints or calibrations.

TTplusWW=TTTo2L2Nu+WWTo2L2Nu+tZq; VZ=WZTo2Q2L+ZZTo2Q2L, preserving the
producer's existing Data-TT-VZ shorthand. This is diagnostic, not deployment.

Reproduce from AlphaRatio_examination with the unbin Python:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python plot_sr_methods.py \
  --pt plots/pt_reweight_2018_channel \
  --mass plots/mass_reweight_2018_channel_fineAlpha \
  --output plots/SR_method_comparison_2018_NEW
```

Supply the original result folders when reproducing from this plot package;
their full source archives are not duplicated here. Hashes record provenance.
