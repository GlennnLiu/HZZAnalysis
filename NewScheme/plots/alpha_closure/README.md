# Data/background alpha-method compatibility, 2018 SR and VR

Compare observed target Data with the full background prediction:

`B_i = alpha_i * (Data_source - TTplusWW_source - VZ_source) + TTplusWW_target + VZ_target`.

Only DY ParticleNet changes, through alpha_i. All other MC retains its identical
full nominal weights. TTplusWW includes TT, WW and tZq; VZ includes WZ and ZZ.
No ZH subtraction or DY pT/mass correction. These are nominal alpha diagnostics,
not deployed production templates. Original alpha ROOT files are untouched.

SR source is full CR, including VR and 105–135. VR source is CR with strict
135 < m2q <165 excluded event by event, so target VR observations are not reused.
Hadronic acceptance is 40 <= m2q <180 and search mass >=300. Target Data is read
and plotted in both SR and VR as requested. Use finalized new reduced trees.

Each ee/mumu channel and six-category selection retains the previous paired-alpha
archive's original channel-specific binning. No channel pooling or bin merging.
The first search-overlapping bin starts at 300; subsequent original edges are
unchanged. Coverage ends at the original alpha histogram's finite last edge.
Yields refer to those covered bins with defined alpha; excluded Data yield from
undefined-alpha bins is explicitly reported, rather than assigning a zero estimate.

## Files

- `six_categories_{ee,mumu}_{SR,VR}_alpha_DataBkg_PNet.{png,pdf}`: four overviews.
- `{category}_{tag}_{channel}_{region}_alpha_DataBkg_PNet.{png,pdf}`: 24 individual plots.
- `all_{SR,VR}_alpha_DataBkg_PNet.pdf`: two twelve-page bundles.
- `compatibility.csv`: 48 channel/method rows plus 24 combined-channel test rows.
- `histograms.npz`: primitives, alpha, prediction variance/covariance and common chi2 masks.
- `metadata.json`: source/archive location and SHA256 plus statistical convention.
- Code and tests retained alongside outputs.

The combined-channel chi2 is the sum of separate ee/mumu chi2s with summed ndof;
its prediction yield sums separately transferred channel estimates. This is NOT
a combined alpha times combined CR yield, and requires no common channel grid.

## Statistical treatment

Write Z = Data_source - other_source and O = other_target. Source and target
are disjoint in each comparison. First-order statistical propagation is:

`Var(B_i) = alpha_i^2 Var(Z) + Z^2 Var(alpha_i) + Var(O)`.

`Cov(B_0,B_1) = alpha_0 alpha_1 Var(Z) + Z^2 Cov(alpha_0,alpha_1) + Var(O)`.

Var(Z) sums the source Data count and other-background sumw2. The paired-alpha
archive supplies correlated DY statistics. The correlated shift variance is
`Var(B_0)+Var(B_1)-2Cov(B_0,B_1)`, not the sum of independent prediction errors.
Second-order products of uncertainty variances are not included. Normalizations
are treated fixed; calibration/modeling/luminosity systematics are not included.

Goodness-of-fit: sum `(Data_target-B_i)^2/(Data_target+Var(B_i))`, with no fitted
normalization. Both choices share a common mask: target Data>=5, DY Neff>=5 in
both alpha regions and both choices, defined alpha, positive total predictions.
ndof is the number of retained bins. Sparse bins remain in plots/yield tables.
Approximate Gaussian chi2 p-values are diagnostics, not complete likelihood tests.
When no bins pass the mask, chi2=0/ndof=0 in CSV is an empty sum, not a fit result;
the readable `compatibility.md` labels this N/A. Merged VR has undefined-alpha
bins excluding 4 untag, 8 btag and 1 VBF target Data events from its yield ratios.
The two predictions are correlated; a difference of their chi2s is not an
independent significance test. SR and VR tests also share some inputs and must
not be summed as independent tests.

Bands show prediction statistical uncertainty. Data markers use Garwood errors;
ratio markers show Data uncertainty only, prediction bands separately around one.
Ratio panels zoom to 0–2.5 with upward markers for off-scale central ratios;
all values are retained in the archive. Signed estimates are not clipped in the
calculations; nonpositive predictions cannot be displayed on log axes.

## Reproduction

```bash
cd /cms/user/liugeliang/CMS/HZZAnalysis/Plot/alpha_closure_pnet_2018
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python -m unittest test_closure -v
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python check_alpha_closure.py \
 --input /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/RedTrees_20261009_pnet \
 --alpha-archive ../alpha_pnet_comparison_2018/histograms.npz --output /path/to/fresh/plots
```

For style-only updates, run `replot_cached.py` using the saved prediction archive.
`validate_package.py` checks figures and propagated statistics and regenerates
the human-readable compatibility summary. Neither needs to reread ROOT inputs.
