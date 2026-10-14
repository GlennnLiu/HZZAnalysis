# SR Data/MC compatibility: DY ParticleNet on/off

2018, six topology/tag categories, ee+mumu combined. Finalized inputs:
`FitTemplate/RedTrees_20261009_pnet`. SR uses strict 70 < m2q <105 GeV and
m2l2q >=300 GeV. Observed SR Data is read and plotted as explicitly requested.

Compare Data with DY+TTplusWW+VZ. Only DY varies between `weight_no_pnet` and
`weight`; non-DY always uses its full nominal weight. TTplusWW contains TT,
WW and tZq; VZ contains WZ and ZZ. No alpha prediction, DY pT/mass reweighting,
ZH subtraction, fitted normalization or new DY DeepJet correction is involved.

`six_categories_SR_m2l2q_DataMC_PNet.pdf` and the analogous `mqq` file show
the overview. Twelve individual figures are saved as PNG/PDF, with a twelve-page
`all_SR_DataMC_PNet.pdf`. Full numerical results: `compatibility.csv` and
`results.json`. Saved bin sums, variances, covariance and chi2 masks:
`histograms.npz`. Histogram row order is noSF yield, withSF yield, noSF variance,
withSF variance, paired covariance (for non-DY the two are identical).

Each goodness-of-fit statistic is:

`chi2 = sum((Data - totalMC)^2 / (Data + sumw2(totalMC)))`.

Use one shared bin mask for both choices: Data>=5, positive total MC, and
effective MC entries>=5 for both. No fitted parameters: ndof is the number of
retained bins. Excluded sparse bins remain in plots and full-range yield tables.
Chi2 includes rate and shape. P-values are Gaussian chi2 approximations, not
full likelihood results. Yield pulls use the complete spectrum with the same
statistical Data+MC variance; Data/MC is also given for each category.

MC bands are statistical only. Data plotting errors use Garwood Poisson intervals;
ratio point errors include Data only, with MC uncertainty shown separately as
bands around one. Calibration/modeling/luminosity/other systematic uncertainties
are not included in chi2 or bands. This is a diagnostic, not a final exclusion
test or a criterion to reject an independently justified calibration. The two
alternatives are highly correlated: their chi2 difference is not an independent
significance test. Likewise mqq and m2l2q use the same events and are not
independent goodness-of-fit tests.

Mass bins match the preceding DY-only comparison: variable m2l2q grid, displayed
as Events /100 GeV. The final bin contains every mass >=3000 GeV and is displayed
3000–4000. Hadronic mass uses 5-GeV SR bins. No high-mass events are cut off.
Ratio panels are zoomed to 0–2.5 to keep the populated spectrum readable;
upward triangles mark central ratios above the range. Sparse-bin uncertainty
stems and bands may extend outside this display range. Full values/errors remain
in the upper panel and archive; no bins are dropped from the yield tables.
Resolved DY has no ParticleNet correction, so its two predictions are identical.

Reproduce in this remote directory with the `unbin` executable:

```bash
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python -m unittest test_sr_pnet -v
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python check_sr_pnet.py \
  --input /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/RedTrees_20261009_pnet \
  --output /path/to/fresh/plots
```
