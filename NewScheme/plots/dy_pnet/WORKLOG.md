# DY ParticleNet comparison worklog — 2026-10-09

Full all-era reduced-tree production completed and aggregate validation passed
for 576 ROOT trees. All comparisons below use those new finalized DY trees, never
partial/pilot outputs or observed Data. The previous plots remain unchanged.

Created `plot_dy_pnet.py`: four eras, six topology/tag categories, ee+mumu combined.
Hadronic mass combines CR and SR only, avoiding double-counting VR already in CR.
m2l2q is plotted in CR, SR and VR separately. Both curves use paired identical
events with full and ParticleNet-excluded weights. Require m2l2q >=300 GeV and
40 <= m2q <180 GeV; full-CR bridge 105–135 is retained.

Statistical histogram variances are sumw2; ratio covariance is the same-event
sum of full-weight times no-ParticleNet-weight. Calibration systematics are not
included. Three regression tests pass, including exact constant-scale cancellation.
Resolved distributions are checked to be exactly unchanged. Signed MC weights
are retained. Variable total-mass bins display events per 100 GeV; the final bin
retains every mass >=3000 GeV, displayed through 4000 GeV without a tail cut.

Saved histograms and yields support restyling without tree reads using
`--histograms histograms.npz`. The final rendering uses density-aware log-axis
limits so a wide tail bin is not hidden by a limit based on unscaled yields.

2018 merged SR integrated yields, combined channels:

| Tag | No ParticleNet | With ParticleNet | Relative change |
| --- | ---: | ---: | ---: |
| untag | 1220.9981 | 1329.6630 | +8.8997% |
| btag | 212.9599 | 177.9716 | -16.4295% |
| vbftag | 116.9209 | 123.2390 | +5.4037% |

CR SF/no-SF ratios are 1.089394, 0.835509 and 1.057843 respectively; SR ratios
are 1.088997, 0.835705 and 1.054037. This shows similar integrated rate changes,
not a proof that the mass-dependent alpha ratio is invariant. No new alpha ROOT
files are computed or deployed here. Shape effects are shown by binwise ratios.

Package: `Plot/dy_pnet_comparison_allEras`, 96 individual PNG/PDF pairs,
16 six-panel PNG/PDF pairs and four 24-page PDF bundles. Includes `yields.csv`,
`histograms.npz`, metadata, code, tests and a complete package validator.
