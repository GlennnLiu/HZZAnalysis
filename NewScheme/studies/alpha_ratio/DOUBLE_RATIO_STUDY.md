# Double ratios and sideband-only transfer corrections

All m(2l2q) >= 300 GeV belongs to the search domain. This study defines no
control region in total mass. A 300--400 GeV reference bin normalizes one
diagnostic; it is not assumed to be signal-free. Corrections are calibrated
only in hadronic-mass sidebands, with both SR and VR excluded from calibration.
SR observations are used for exploratory validation, never to determine the
central correction. A signal in SR could affect that diagnostic. Potential
signal leakage into hadronic sidebands requires explicit treatment in a search.

## Reproduce

```bash
source /scratchfs2/cms/liugeliang/miniconda3/bin/activate
conda activate unbin
cd /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate
python AlphaRatio_examination/double_ratio_untag.py \
  --base RedTrees --cr-base RedTrees_old \
  --output AlphaRatio_examination/plots/double_ratio_untag_2018 --bootstrap 500
python AlphaRatio_examination/plot_double_ratio_summary.py \
  AlphaRatio_examination/plots/double_ratio_untag_2018
```

Inputs are broad CR trees plus current SR trees; CR already contains VR.
Therefore separate VR trees are not concatenated. Nominal transfer bins are
300,400,500,650,850,1100,1600,3500 GeV in resolved; merged uses
300,500,850,3500 GeV because of lower statistics. Hadronic bins are
40,55,70,85,105,120,135,150,165,180 GeV. The data and MC weights are inherited
from existing reduced trees. Other-background systematic uncertainties and
event-level production weight correctness are not established by this study.

## Double ratios

Let Z(x,q)=Data-TTplusWW-VZ and D(x,q)=DY MC; r=Z/D.
The four-region diagnostic is

```
DR(x,q) = r(x,q)/r(x_reference,q) / [r(x,q_reference)/r(x_reference,q_reference)]
```

The references are the first total-mass bin and 40--55 GeV hadronic mass.
If r factorizes into f(x)*g(q), this quantity equals one. The reference row and
column equal one algebraically; their zero errors do not imply precise
measurements. Shared-reference statistical covariance is included in the
global tests. Tests exclude cells with effective statistics below four;
high-mass heatmap colors must be interpreted with their printed errors.

There is also a reference-free-in-total-mass diagnostic:

```
DR_SR_VR(x) = [Z_SR(x)/Z_VR(x)] / [D_SR(x)/D_VR(x)]
```

This is the ratio of SR and VR alpha closure factors, since their common CR
cancels. It checks directly whether a VR-derived correction could be applied
unchanged to SR. It is a diagnostic involving SR data, not a calibration rule.

Combined resolved four-region tests (statistical errors only):

| Cells tested | chi-square / rank | p-value |
|---|---:|---:|
| hadronic sidebands | 16.23 / 20 | 0.702 |
| SR | 22.35 / 10 | 0.0134 |
| VR | 6.11 / 10 | 0.806 |

Thus the sidebands and VR do not demand strong nonfactorization at this
precision. SR gives moderate tension; this does not prove the origin is DY.
The tests use Gaussian/log delta-method approximations and are not search
significances. No look-elsewhere adjustment or systematic covariance is applied.

## Candidate correction

Fit log r(x,q) in 40--70, 105--135, and 165--180 GeV hadronic sidebands only.
The simplest model is

```
log k(x,q) = a0 + a1*u + a2*u^2 + b1*v + b2*v^2
u=(q-100)/60; v=log(x/500)
```

This factorizes into separate total-mass and hadronic-mass corrections. Other
tested candidates add u*v, u^2*v, or u^3, or interpolate locally between five
sideband slices after smoothing their total-mass dependence. Central predictions
use bin-center approximations to the fitted surface and DY-weighted sums over
hadronic bins; a production implementation should integrate over actual event
kinematics and verify binning sensitivity.

For target T and sideband union C, the proposed alpha correction is

```
F_T(x) = <k(x,q)>_DY,T / <k(x,q)>_DY,C
alpha_corrected_T(x) = alpha_nominal_T(x) * F_T(x)
prediction_T(x) = alpha_corrected_T(x) * Z_C(x)
```

No observed Z_SR or Z_VR enters F. A total-mass-only correction cancels in F;
the hadronic correction changes F because DY's q distribution differs by region
and total mass. The empirical correction need not be identical in SR and VR.

Combined resolved held-out validation, in seven total-mass bins:

| Model | SR observed/prediction | SR statistical chi-square | VR observed/prediction | VR statistical chi-square |
|---|---:|---:|---:|---:|
| nominal broad-CR alpha | 0.9752 | 66.24 | 0.9171 | 167.45 |
| separable quadratic | 0.9721 | 26.79 | 1.0114 | 6.16 |
| linear interaction | 0.9721 | 25.97 | 1.0114 | 5.57 |
| quadratic interaction | 0.9736 | 20.28 | 1.0101 | 5.38 |
| cubic separable | 0.9488 | 36.65 | 1.0319 | 13.45 |
| local interpolation | 0.9566 | 38.72 | 0.9983 | 2.52 |

Prediction statistical covariance is obtained from 500 Gaussian cell-yield
resamples, refitting the calibration on each resample. This preserves dependence
between the sideband fit and the sideband yield to which alpha is applied.
The resampling approximates weighted-yield statistics, not exact event bootstrap.
Negative signed cells are not clipped; log fits omit unstable cells. Sparse
merged fits, especially local interpolation, lose resamples and should not be
used for a precision correction. Covariance includes finite MC statistics;
background-subtraction and modeling systematics are absent. Covariance is model
conditional and does not include the interpolation-model choice.

Cross-channel test with the simplest correction:

| Sideband training | Validation channel | SR closure | VR closure |
|---|---|---:|---:|
| ee | mumu | 0.9694 | 1.0150 |
| mumu | ee | 0.9753 | 1.0099 |

This supports consistency across lepton channels, but neither lepton channel
is a signal-free control sample. Both tests hold target-region data out of
calibration. Separate-channel fits have independent event statistics but shared
physics/systematic assumptions.

Resolved subtraction stress tests, refitting sidebands after each variation:

| TT/WW+tZq scale | VZ scale | SR closure | VR closure |
|---:|---:|---:|---:|
| 1.0 | 1.0 | 0.97207 | 1.01138 |
| 0.8 | 1.0 | 0.97214 | 1.01127 |
| 1.2 | 1.0 | 0.97201 | 1.01149 |
| 1.0 | 0.7 | 0.97974 | 1.01059 |
| 1.0 | 1.3 | 0.96431 | 1.01218 |

The 20% and 30% variations are exploratory stress tests, not established nuisance
uncertainties. Common TT normalization changes largely cancel here. VZ affects
SR more because it contains a hadronic-Z contribution; a VZ normalization/shape
audit is warranted before claiming the remaining mismatch is solely DY.

## Recommendation for the search

Retain the broad untagged hadronic CR. Use the simplest sideband correction
as a candidate for resolved categories: it strongly improves VR closure and
reduces the SR shape mismatch without calibrating on SR. It is not a complete
solution: about 3% integrated SR nonclosure remains, with a mass-dependent
residual. Prefer a constrained transfer-factor model and an explicit SR
nonclosure nuisance to an exact point-by-point correction of observed SR data.
The size and correlation of that nuisance need independent modeling checks,
closure studies, and signal-injection tests; the observed SR residual alone is
not an adequate search calibration or a justified universal 3% uncertainty.

For merged categories, keep nominal alpha for now: nominal SR closure is
already satisfactory and the sideband corrections worsen central yield closure.

The final likelihood should include signal in every fitted mass slice. Schematically,

```
expectation_SR = alpha_corrected * B_DY_CR + B_TT_SR + B_VZ_SR + mu*S_SR
expectation_CR = B_DY_CR + B_TT_CR + B_VZ_CR + mu*S_CR
```

Use CR observations once through this likelihood, rather than fixing a
data-derived prediction and treating its input CR as an independent measurement.
The q-sideband correction calibration also reuses CR events: joint fitting or
retaining their covariance is necessary. If signal leakage is ignored, the
effective reconstructed signal after subtraction can be S_SR-alpha*S_CR.
No total-mass sideband is assumed or available in this wide-mass search.

## Artifacts

`plots/double_ratio_untag_2018` contains PNG/PDF heatmaps, SR/VR double-ratio
plots, all-model correction plots, and clean `*_candidate_summary` plots.
`cells.csv`, `double_ratios.csv`, `region_double_ratios.csv`,
`double_ratio_tests.csv`, `validation.csv`, `fit_quality.csv`, and
`correction_factors.csv` retain numerical results and candidate factors.
The scripts and plots are exploratory; production alpha files are unchanged.
