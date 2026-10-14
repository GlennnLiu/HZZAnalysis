# 2018 SR Data/MC compatibility — 2026-10-09

Explicit user request to inspect observed SR Data. Source is the completed,
validated `RedTrees_20261009_pnet` production. No further tree production,
alpha deployment or modifications to existing plotting outputs.

Six categories, ee+mumu combined; m2l2q >=300 GeV and strict SR m2q 70–105.
Full prediction is DY+TTplusWW+VZ. Only DY switches ParticleNet on/off;
all other backgrounds retain the same nominal full weights. No ZH subtraction.
Added mass overlays, Garwood Data errors, MC statistical bands and Data/MC panels.
Ratio range is 0–2.5 with triangles for off-scale central points; sparse tail
values/uncertainties are retained in archives and full-spectrum yield summaries.

Chi2 uses Data+MC statistical variance and identical bin masks for both choices:
Data count >=5, positive total MC, MC Neff >=5 for both choices. No fitted
normalization. All excluded sparse bins remain plotted and included in yields.
P-values are Gaussian approximations, with no modeling/calibration systematics.
Three unit tests passed for chi2, masks and zero-count Poisson upper intervals.

| Category | Data/MC noSF | Data/MC withSF | m2l2q chi2/ndof noSF | withSF |
| --- | ---: | ---: | ---: | ---: |
| Resolved untag | 0.9448 | 0.9448 | 252.60/12 | 252.60/12 |
| Resolved btag | 0.9366 | 0.9366 | 12.14/7 | 12.14/7 |
| Resolved VBF | 1.0140 | 1.0140 | 24.54/10 | 24.54/10 |
| Merged untag | 1.1471 | 1.0791 | 37.72/11 | 19.91/11 |
| Merged btag | 1.3472 | 1.5415 | 26.14/9 | 45.44/9 |
| Merged VBF | 1.1185 | 1.0741 | 9.26/8 | 8.81/8 |

Hadronic-mass chi2/7 for merged untag improves 31.82 ->11.94; btag worsens
26.77 ->45.74; VBF improves 3.36 ->2.49. These tests share events with m2l2q
and are not independent. SF improves merged untag, slightly helps VBF and
worsens btag: not a uniform remedy for direct-DY nonclosure. An independently
justified calibration should not be accepted/rejected merely by this closure
diagnostic. No significance is assigned to differences between correlated chi2s.

Retained code, tests, input-path/method documentation, 24 numerical rows,
bin archives, 12 individual PNG/PDF pairs, two overview PNG/PDF pairs and a
twelve-page PDF bundle. Cached rerendering changes display only, not statistics.
