# DY with/without ParticleNet SF, all eras

Source: validated `FitTemplate/RedTrees_20261009_pnet` production. The full
576-tree production passed aggregate validation on 2026-10-09.

Remote package: `/cms/user/liugeliang/CMS/HZZAnalysis/Plot/dy_pnet_comparison_allEras`.
Plots combine ee and mumu, separately for each of the six topology/tag categories
and each era. Both curves use precisely the same events, comparing `weight` with
`weight_no_pnet`. Other corrections, signed MC weights and normalization are
unchanged. ParticleNet is applied only to merged events; resolved overlays and
ratios are validated to be identical and one, respectively.

## Files

- `{era}_{category}_{tag}_mqq`: hadronic mass in CR + SR, without double-counting
  VR (which is already part of full CR).
- `{era}_{category}_{tag}_{CR,SR,VR}`: m2l2q in the named hadronic-mass region.
- `six_categories_{era}_{mqq,CR,SR,VR}_DY_ParticleNet`: six-panel overview pages.
- `all_{era}_DY_ParticleNet.pdf`: 24 individual plots per era, bundled.
- `yields.csv`: yields, integrated SF/no-SF ratios and correlated statistical errors.
- `histograms.npz`: unscaled yields, variances and same-event covariance in every bin.
- `metadata.json`: input path, selection, statistical treatment and tail convention.
- `plot_dy_pnet.py`, `test_dy_pnet.py`: retained code and three passing ratio tests.

Each figure is saved separately in PNG and PDF. Absolute-yield comparisons,
not unit-normalized shapes: changes in both rate and shape are visible.

## Selection and uncertainty conventions

Both masses require m2l2q >= 300 GeV and 40 <= m2q < 180 GeV. CR retains
105–135 and VR; SR uses the producer's strict 70 < m2q < 105 bounds, and VR
strict 135 < m2q < 165. No observed Data is read, including in SR.

Hadronic mass has 5-GeV bins. m2l2q uses a common explicitly saved variable-bin
grid, shown as Events / 100 GeV. All masses >= 3000 enter the final displayed
3000–4000 bin, with density based on its 1000-GeV display width. This is not a
cut at 4000 GeV; yield tables retain the complete tail.

Shaded bands are MC statistical uncertainties, sqrt(sum(weight^2)), not ParticleNet
calibration systematics. For ratio R = A/B, with A the SF-weighted histogram,
the variance includes shared-event covariance C = sum(weight*weight_no_pnet):

`Var(R) = (Var(A) + R^2 Var(B) - 2 R C) / B^2`.

Ratios are omitted when the signed denominator is nonpositive; their histogram
contents are not clipped or discarded in the archive. A constant event SF has
zero ratio statistical error because the two histograms are fully correlated.
These error bars measure uncertainty from the finite sampled SF/flavour mixture,
not uncertainty of the SF calibration itself. CR and VR plots are overlapping
selections and must not be treated as independent observations.

## Reproduction

```bash
cd /cms/user/liugeliang/CMS/HZZAnalysis/Plot/dy_pnet_comparison_allEras
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python -m unittest test_dy_pnet -v
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python plot_dy_pnet.py \
  --input /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/RedTrees_20261009_pnet \
  --output /path/to/fresh/plot/directory
```

No production alpha ratios, reweighting fits or old plots are modified.

For a style-only replot, add `--histograms histograms.npz` to the plotting command;
this uses the saved bin sums/covariance and avoids all ROOT-tree reads. Run
`validate_package.py` from the retained package after generation to check every
figure, histogram and exact resolved-category equality.
