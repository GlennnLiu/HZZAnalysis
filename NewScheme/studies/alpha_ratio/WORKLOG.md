# Alpha-ratio investigation worklog

Updated 2026-10-06, Europe/Paris. Remote project:
`/cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/AlphaRatio_examination`.
This directory is mirrored locally under `LXLOGIN/CMS/HZZAnalysis/FitTemplate`.

## Scope and earlier investigation

The search covers all four-body masses >=300 GeV. There is no signal-free
control region in total mass. CR/SR/VR in these studies describe hadronic mass.
The user prefers retaining 105--135 GeV in the broad hadronic CR and subtracting
ZHbb, with category-dependent validation.

Existing detailed records are preserved:

- `CLOSURE_STUDY.md`: untagged broad/narrow CR closure. Broad CR improves the
  resolved closure, but residual mass dependence remains; alpha cancels common
  DY normalization mismodeling, not region-dependent mismodeling.
- `DOUBLE_RATIO_STUDY.md`: sideband-only transfer-model candidates. Useful
  diagnostics, not a deployed production correction or a total-mass CR.
- `AUXILIARY_ALPHA_STUDY.md`: the user's intended independent-selection double
  ratio. MET>150 is small and top dominated; resolved weights outside nominal
  DeepJet validity are missing. Cached auxiliary unit-SF results are provisional.
  No precise or transferable correction has been established. Flipping
  foundZ1LCandidate while retaining ee/mumu produces an empty sample.

## Completed ZHbb skim and subtraction study

The isolated 2018 ANATree+DeepJet workflow is complete. See
`../../ANATree/zhbb_2018/WORKLOG.md` and its README/validation files. Final skim:
`ANATree/zhbb_2018/ZH_HToBB_ZToLL.root`, with 496,197 retained entries from
4,837,847 processed, normalization 4,547,627 and 356 scalar branches unchanged.
The analysis card's 0.0317 pb and 59.83 fb^-1 are inherited conventions.

`prepare_zh_reduced.py` created 36 isolated ZH reduced trees under
`zhbb_reduced_2018`. `compare_zh_subtraction.py` compared all six standard
categories, ee/mumu/combined, SR/VR, using broad-CR `RedTrees_old` and current
target `RedTrees`. The broad CR excludes both SR and VR and retains 105--135.
Nominal DY ratios were recomputed in common coarse bins, including the full
infinite tail, with deterministic DY-only merging of nonpositive DY bins.
This is not an exact reproduction of stored fine-bin production alpha.

Background definitions and fair comparisons:

```
B = TTplusWW + VZ
TTplusWW = TTTo2L2Nu + WWTo2L2Nu + tZq
VZ = WZTo2Q2L + ZZTo2Q2L
P0 = alpha * (Data_CR - B_CR)
P1 = P0 - alpha * ZHbb_CR
O0 = Data_target - B_target
O1 = O0 - ZHbb_target
```

Compare P0/O0 and P1/O1 for DY closure. To use the same Data-TT-VZ target
points for both predictions, compare P0 and P1+ZHbb_target against O0.
Do not confuse the DY-only shift with the complete-background shift.
Before/after predictions are correlated; common baseline fluctuations cancel
in the difference covariance.

The study passed its analytic regression tests and output validation, with a
finer-bin cross-check. Nominal DeepJet applies only to resolved non-VBF btag/untag.
Merged Hbb uses the existing ParticleNet convention, not a newly validated Hbb
calibration; unit-tag checks are saved. Only ZH_HToBB_ZToLL is included, not
ZHcc, ggZH or every Higgs process. Production alpha/configuration/subtraction
inputs have not been updated to consume these isolated reduced trees.

## Statistical coverage result

Combined ee+mumu integrated DY shifts compared with the ORIGINAL prediction
statistical uncertainty (all masses >=300, full tail):

| Category | SR | VR |
| --- | ---: | ---: |
| resolved btag | 0.313 sigma | 0.331 sigma |
| resolved untag | 0.044 sigma | 0.041 sigma |
| resolved vbftag | 0.012 sigma | 0.012 sigma |
| merged btag | 0.139 sigma | 0.153 sigma |
| merged untag | 0.019 sigma | 0.016 sigma |
| merged vbftag | 0.009 sigma | 0.012 sigma |

In the primary binning, all 246 SR/VR category/channel DY-only bin shifts lie
inside the original statistical band; the largest is about 0.824 sigma in
merged-btag ee VR >=1600. SR alone has maximum 0.716693 sigma in merged-btag ee
>=1600; combined SR maximum is 0.600175 sigma in that category/tail.
The finer-bin cross-check also stays below one sigma for the DY-only shift.

These statements concern P1-P0 relative to sqrt(Var(P0)), not the significance
of the difference, whose shared-error variance is smaller. Statistical errors
include signed CR sumw2 and finite DY transfer statistics, not tagging,
luminosity, normalization or modeling systematics. A known bias can and should
be subtracted even when below current statistical precision.

The complete-background shift P1+ZHbb_target-P0 has a sparse VR exception
when judged only against P0's band: merged-btag ee >=1600, with zero CR data
and a small negative signed P0. It is not an SR DY-only exception and remains
small compared with the target-data comparison uncertainty. Do not generalize
the DY-only coverage statement to every definition of the plotted prediction.

ZH subtraction cannot account for the resolved-untag residual: it changes the
SR complete prediction by only about +0.885 events against a roughly 1,292-event
baseline excess, and VR by about -2.996 against roughly 3,825. No unique cause
of the residual or production correction has been established.

## SR-only uncertainty-band figures and CMS restyling

The user requested P0 with its +/-1 sigma statistical band, P1 on the same
plot, and no SR Data-TT-VZ observations. `plot_sr_zh_statband.py` reads only
edges, pred_nozh, pred_withzh and transferred_zh from the saved archive;
target-observation arrays are not read. P1 is DY-only: target ZH is not added.
The lower panel shows (P1-P0)/sqrt(Var(P0)) with a [-1,+1] band.

The initial plots are preserved in `plots/zh_subtraction_2018_SR_statband/`.
The current presentation is `plots/zh_subtraction_2018_SR_statband_CMS/`:

- `combined_SR_ZH_statband_summary.pdf` / `.png`: all six combined-channel categories.
- `<resolved|merged>_<btag|untag|vbftag>_<ee|mumu|combined>_SR_ZH_statband.pdf`
  and `.png`: 18 separate category/channel figures.
- `all_SR_ZH_statband.pdf`: 19-page bundle, overview first.
- `summary.csv`, `metadata.json`, `README.md`, and a reproduction-script copy.
- Sibling `.zip`: all separate files and the bundle. Files are on LXLOGIN and local.

HZZAnalysis already uses mplhep CMS style in make_1Dtemplate.py; ROOT CMS_lumi
helpers also exist in SigRes/SigEff/Templates2D. The installed mplhep is 0.3.23.
Restyled figures use `hep.style.CMS`, inward ticks, CMS Work in progress and
59.8 fb^-1 (13 TeV), rounded only for display. They are not labelled Simulation
because the estimates use CR data. They do not assert CMS public approval.

Removed the original high-mass hatch: it only denoted an open-ended bin, not
exclusion or blinding. The final bin is now labelled >= its lower edge, with an
overflow note. Only statistical bands are shaded; signed values and all tail
events are retained. Each figure is separately saved in PNG and PDF.

Verification: restyled summary.csv is byte-identical to the original; source
archive digest is identical; 19 PNGs and 20 PDFs (including bundle) are present.
The six-category overview and representative merged-btag PNG were visually
checked for layout. Restyling changed no yields, variances, alpha or selections.

## Saved artifacts and reproduction

Full SR/VR closure study: `plots/zh_subtraction_2018_broadCR/`, 66 figures each
in PDF/PNG, 36 yield rows, bins.csv, histograms.npz, metadata, audits and binning.
Finer-bin cross-check: `plots/zh_subtraction_2018_finecheck/` (also included in
the broad-study archive). Source histogram SHA256:
`d7314a9ab0b7b611645add21f4f76636dc01750e767c0b1cd48cf63439df111a`.

Archive arrays pred_nozh/pred_withzh/transferred_zh have shape (2,n): row 0 is
sumw, row 1 variance. Keys use category__tag__channel__target__quantity.
Combined rows sum the separately transferred ee/mumu predictions and variances.
Edges end at +infinity; never drop the last bin during plotting or integration.

```bash
cd /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/AlphaRatio_examination
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python plot_sr_zh_statband.py \
  --input plots/zh_subtraction_2018_broadCR/histograms.npz \
  --output plots/zh_subtraction_2018_SR_statband_CMS_NEW
```

Use a fresh output path: the plotter refuses existing nonempty directories.
No skimming, ROOT import or batch submission is needed for restyling.

## Handoff and remaining work

Completed: skim, isolated reduced-tree production, category-dependent SR/VR
comparison, statistical-coverage assessment and CMS-style SR-only plots.
Not completed or authorized by these diagnostics: deployment in production,
a validated auxiliary double-ratio correction, independent Hbb tag calibration,
or a definitive explanation of the remaining DY closure discrepancy.

Documentation-only update on 2026-10-06: added this worklog, the skim worklog,
and a focused alpha-ratio reference/routing update to cms-hzz-fittemplate.
Preserved old study records, old figures, analysis files and unrelated changes.

## 2026-10-07: DY hadronic-pT reweighting

Completed the user's new test with no ZH subtraction: fit
(Data-TTplusWW-VZ)/DY versus pt2jet (resolved) or ptmerged (merged) in C+VR,
then apply those factors to DY SR and recompute nominal/reweighted alpha.
Both category-specific and topology-shared tag corrections are tested; ee/mumu
share each fit and are compared separately and combined. CR-only SF fits provide
VR-held-out validation. For SR the primary alpha source is C+VR; for VR it is
C alone, so the target never overlaps its own alpha source. SR C-only alpha is
also saved for comparison with the previous study's convention.

The existing reduced trees lack hadronic-pT branches. `pt_reweight_scan.py`
therefore used the current `make_reduced_trees.py` input paths, selections and
MC weighting helpers on original combined skims, without changing those trees
or the producer. One chunked scan of all 11 enabled samples is cached in
`pt_reweight_source_2018`. All 144 process/category/channel/region sumw/sumw2
comparisons match the previous reduced-tree study (maximum relative difference
7.25e-14). All search masses >=300, including the full tail, remain in scope.

`pt_reweight_core.py` propagates analytic first-order Jacobians through the SF,
weighted DY and alpha construction, including shared calibration/CR/channel
correlations and VR reuse. Seven regression tests passed, including numerical
gradient checks and rejection of an alpha source containing its target.
No SR observations enter the SF fit. Statistical-only Gaussian diagnostics
exclude inadequately populated bins from chi2 but retain them in distributions
and tables. Tag/subtraction/modeling systematics and adaptive-bin-choice errors
are not established by this test.

Final outputs: `plots/pt_reweight_2018_fullCR/`, 89 separate CMS-style PDF/PNG
pairs, a 71-page search-mass/SF/held-out-VR bundle, an 18-page SR-pT bundle,
numerical results, covariances, source cache/audit, validation and reproduction
scripts. The README is the detailed report. The sibling ZIP is the full package;
both package and folder are mirrored locally and on LXLOGIN. An intermediate
`plots/pt_reweight_2018` used C-only SR alpha and is not the final presentation.

Key category-specific ee+mumu untagged results:

- Resolved direct SR closure: 0.9409 ->0.9951, but search-mass chi2/rank remains
  119.85/10 after reweighting (253.38/10 before). Normalization alone gives
  0.9965 and 127.17/10: much of the improvement is not a pT-shape effect.
- Resolved SR alpha with C+VR source: 0.9995 ->0.9982, 65.42/10 ->62.54/10.
  Good integrated yield closure does not imply mass-shape closure.
- Resolved VR alpha with C source: 0.9170 ->0.9192, 167.24/10 ->159.66/10.
- Merged direct SR: 1.2079 ->0.9786, 29.80/4 ->1.34/4; normalization alone
  already gives 0.9750 and 1.60/4. Its SR alpha changes 0.9750 ->0.9784,
  while VR alpha changes 1.2282 ->1.2126; neither changes appreciably.
- With the old C-only SR-alpha convention, resolved untagged changes
  0.9753 ->0.9746. Do not credit the switch to a C+VR source as a pT benefit.

The identity alpha_rw/alpha_nominal = <r>_DY,target / <r>_DY,source explains
why a common normalization cancels and pT corrections largely cancel when
the conditional source/target pT distributions are similar. This test improves
direct DY normalization but does not resolve the transfer-model discrepancy.

Resolved category-sharing calibration has chi2/rank 31.81/17 (p=0.0159),
with the VBF signed calibration yield 1.051 +/-0.013 relative to a common SF.
Shared SF worsens resolved VBF direct SR closure (1.015 ->1.069), whereas its
own SF gives 1.017. Do not automatically pool resolved VBF with untagged.
Merged pooling is compatible at the current statistical precision (13.87/13,
p=0.383), but separate-category curves and sparse-tag uncertainties are retained.

The primary VR SF includes VR and is consequently partly in-sample. C-only
fits still leave the resolved-untagged VR residual, rather than demonstrating a
transfer cure. No pT correction has been deployed into production templates,
and no new batch or skimming jobs were submitted for this study.

## 2026-10-07: calibration-region pT before/after plots

Added `plot_cr_pt_reweight.py` and `CR_PT_PLOTS.md`, reusing the unchanged
`plots/pt_reweight_2018_fullCR/source_histograms.npz` cache and original
adaptive pT factors. Outputs are isolated in `plots/pt_reweight_2018_CR_closure/`
with a sibling ZIP: 48 separate CMS-style PNG/PDF pairs and a 48-page
`all_CR_pt_reweight.pdf`, plus yield/bin CSVs, covariance NPZ and validation.
Both C+VR primary calibration and C-only cross-check are shown for resolved/
merged, each tag and their pooled sum, ee+mumu and individual channels.

All 16 combined fitted samples pass exact bin-value and Jacobian closure
assertions. For a category-specific factor the combined category closes; for
a pooled factor only the combined pooled sample is forced to close. Separate
channels and tags under a shared fit are not independently forced to unity.
The covariance accounts for reuse of the calibration data: the fitted-sample
ratio error cancels rather than double-counting data and fitted-DY errors.
This is an in-sample construction check, not independent SR validation.
No ZH subtraction, correction change, source scan, ROOT rewrite or batch job.

## 2026-10-07: independent ee/mumu pT corrections per category

At the user's request, added isolated `pt_reweight_channel.py`,
`plot_pt_reweight_channel.py`, `test_pt_reweight_channel.py` and
`PT_CHANNEL_STUDY.md`. Category-inclusive/shared-tag fits and curves are
omitted. Resolved/merged x btag/untag/vbftag x ee/mumu are calibrated separately
in C+VR (primary) and C only (held-out-VR cross-check). Combined-channel
plots sum independently corrected/transferred channel Estimates; they do not
use a combined calibration. No source rescan, skimming or production update.

All 13 tests passed (seven original plus six new). Runtime assertions verify
36 calibration spectra close in values and Jacobians, 36 SR pT/mass projections
agree in integrated values/derivatives, and each factor uses only its own
category/channel with no SR data. Nominal yields/errors match the original
study. A first numerical attempt caught mass-grid sensitivity of nominal alpha;
the final calculation retains the previous mass edges in all six categories,
so changes are attributable to the channel split, not rebinning. No additional
merging was necessary. Sparse pT bins use the original per-fit Neff>=25 rule;
no category or channel pooling. Independent merged-VBF fits retain only two
or three pT bins, reflecting their limited statistical precision.

New output: `plots/pt_reweight_2018_channel/` and sibling ZIP, mirrored locally
and on LXLOGIN. 147 separate CMS-style PDF/PNG pairs: 93-page mass/SF bundle,
18-page SR pT bundle, 36-page calibration pT bundle. Main six-category overviews
exist for ee, mumu and their sum. Every category/channel has nominal/reweighted
direct SR and alpha SR/VR plots plus held-out-VR alpha validation. Calibration
and SF plots show only category/channel-specific corrections. Numerical files,
full covariance, unchanged source cache/audit, prior baseline summary/mass edges
and reproduction scripts are included. No earlier plot folder is overwritten.

Primary untagged channel-split results (observed/predicted, ee+mumu):
resolved direct SR 0.99517 (chi2/rank 119.61/10); SR alpha 0.99819 (62.54/10);
VR alpha 0.91917 (159.74/10). Merged direct SR 0.97811 (1.29/4), SR alpha
0.97844 (1.21/4), VR alpha 1.21209 (6.47/4). Splitting the calibration by
lepton channel barely changes the combined untagged conclusion. Resolved
channel-specific VR alpha closure is 0.92239+/-0.00976 (ee) and
0.91676+/-0.00853 (mumu), statistical only. Primary VR calibration still
includes VR; use the held-out-VR plots for the independent cross-check.
No ZH subtraction. Systematic/modeling uncertainties remain outside this study.

## 2026-10-07: mass distribution in the pT-calibration CR

Added isolated `plot_cr_mass_pt_reweight_channel.py` and `CR_MASS_PT_CHANNEL.md`.
Output `plots/pt_reweight_2018_channel_CR_mass/` plus sibling ZIP: 42 CMS-style
PNG/PDF pairs and a 42-page bundle. Includes six-category overviews and all
ee/mumu/combined distributions in primary C+VR and C-only cross-check. Curves
show nominal DY, channel-specific pT reweighting and per-channel normalization
only, versus Data-TTplusWW-VZ in that calibration region. All existing factors
and mass grids remain unchanged; no ROOT/source/production changes or ZH.

All 13 existing/channel regression tests passed. Recomputed 24 factors match
the saved edges, values AND covariances. All 36 calibration mass spectra have
exact integrated corrected/norm-only yield and Jacobian closure. This does not
enforce bin-by-bin mass closure. Joint residual covariance retains calibration
data reuse and its normalization constraint; rank is derived from eigenvalues.
Numerical tables, full covariances, source cache/audit, reference calibration
and reproduction code are packaged. Original outputs are not overwritten.

Untagged C+VR, combined ee+mumu mass chi2/rank:
resolved nominal 384.62/10, normalization-only 37.58/9, pT-reweighted 70.78/9;
merged nominal 79.84/4, normalization-only 1.45/3, pT-reweighted 0.62/3.
Raw normalization mismatch is reduced, but resolved shape remains statistically
inconsistent and is WORSE with pT correction than with normalization alone.
This holds separately: ee norm 36.42/9 vs pT 57.51/9, mumu norm 11.48/9 vs
pT 23.97/9. Merged passes at this precision/coarse binning after either option.
Integrated yield closure is forced and should not be used as shape validation.
This is an in-sample statistical-only diagnostic, not independent evidence
that the correction transfers to SR or a unique diagnosis of DY mismodeling.

## 2026-10-07: replace pT calibration by CR m2l2q reweighting

Implemented `mass_reweight_dy.py`, `plot_mass_reweight.py`,
`test_mass_reweight.py`, `MASS_REWEIGHT_STUDY.md`. Inputs are the original
nominal DY/signed-Z cells, NOT the pT-reweighted templates. Each resolved/
merged x btag/untag/vbftag x ee/mumu is independently calibrated in C+VR
and C-only; no ZH subtraction or category-inclusive calibration. Source/input
selections, weights, ROOT data and production alpha remain unchanged.

Primary output `plots/mass_reweight_2018_channel_fineAlpha/` plus ZIP, local
and LXLOGIN mirrors. 135 individual CMS-style PNG/PDF pairs, 93-page mass
comparison/SF and 42-page calibration bundles. Includes direct SR, alpha
SR/VR, held-out-VR alpha checks, SFs and CR mass comparisons for all channels.
Prior `plots/mass_reweight_2018_channel` is a completed coarse-alpha-grid
diagnostic, NOT the primary before/after comparison with previous studies.

Correction bins start from old grids and merge until Z/DY Neff>=25 per
calibration/channel and target MC positivity. Number of correction bins:
resolved btag/untag/VBF 4/9/7, merged 2/3/1. No SR data sets the grid or fit.
Primary VR enters calibration; C-only cross-check holds VR out. Sparse merged
VBF can only support a normalization correction. Common boundaries do not
pool any channel or category correction values.

Important refinement: original alpha/comparison grids are fully preserved.
Coarse mass factors are expanded as constants onto constituent alpha bins;
no alpha bin crosses a factor boundary. All nominal SR/VR yields and errors
match the previous channel-dependent study. This avoids confusing a change
in nominal alpha from rebinning with a correction benefit. Fine CR plots need
not close within a broad correction bin; fitted-grid and integrated closures
are separately checked.

15 regression tests passed, including finite differences, isolation, no SR-data
grid choice, nonconstant mass-factor cancellation and coarse-factor expansion.
Runtime verifies 36 fitted-grid calibration value/Jacobian closures, 72 alpha
SR/VR cancellations and 36 matched-source direct/alpha identities on the
correction grid. Full first-order statistical covariance retains calibration
data reuse; alpha cancellation includes fitted-factor uncertainty. Output
covariance matrices are also checked before packaging; tests.log is saved.

Main untagged ee+mumu results on preserved original grids:
resolved direct SR observed/predicted 0.94088 ->0.99951 (+/-0.00633 after),
chi2/rank 253.38/10 ->65.43/10. Merged direct SR 1.20793 ->0.97482
(+/-0.03913 after), 29.80/4 ->1.33/4. Normalization-only gives 127.06/10
and 1.56/4 respectively. Alpha SR remains 0.99951,65.42/10 (resolved)
and 0.97498,1.20/4 (merged); alpha VR remains 0.91705,167.24/10 and
1.22824,6.59/4. Central values AND uncertainties of alpha are unchanged.

Algebra: alpha_rw=(r*DY_target)/(r*DY_source)=alpha_nominal in every original
bin. Direct r*DY_target equals ordinary alpha*Z_source when both use the
correction grid and source. Finer displayed direct templates can differ from
fine-bin alpha because sparse correction bins explicitly smooth the CR ratio.
This improves absolute DY but is not an additional cure for alpha nonclosure.
Statistics-only diagnostics, not production deployment or a final likelihood.

## 2026-10-07: five-method SR overlays and skill handoff

Added plot_sr_methods.py and SR_METHOD_COMPARISON.md in remote analysis root,
with local development copies under method_comparison_2018/. Outputs:
plots/SR_method_comparison_2018/ plus ZIP, 18 individual category/channel
PNG/PDF pairs, 3 six-category ee/mumu/combined overviews and a 21-page bundle.
Each SR plot compares nominal direct DY, direct pT-reweighted DY, direct
mass-reweighted DY, nominal alpha and pT-reweighted alpha with target signed
Data-TTplusWW-VZ. No ZH subtraction, no new fit, no tree/production change.

Use latest pt_reweight_2018_channel and mass_reweight_2018_channel_fineAlpha.
Validation matched primitive-source SHA256 plus 8 fields (original edges,
observation/value covariance, comparison mask, nominal DY/value covariance,
nominal alpha/value covariance) in all 18 category/channel combinations.
All 90 prediction rows, full covariances and joint closure-bin errors are
copied unchanged into results.npz/summary.csv. CMS style, prediction bands,
ratio-panel points and integrated observed/predicted ratios in individual
legends. Combined overview and resolved-untag combined visually inspected.
Near-overlapping methods retained; no independent significance is assigned
to their differences because they share source data and MC.

Updated cms-hzz-fittemplate skill routing/reference with latest correction
artifacts, fineAlpha versus coarse-alpha distinction, input schema/provenance,
calibration definitions, same-bin factor cancellation and covariance caveats.
Skill draft validated with skill-creator quick_validate; plot script compiled.
Local/remote worklogs and local skill reference are maintained. No unrelated
git edits, ROOT inputs, producer settings or previous plots were modified.

## 2026-10-07: five-method VR comparison including direct DY

Added plot_vr_methods.py and VR_METHOD_COMPARISON.md in analysis root,
local copies under method_comparison_2018/. Outputs plots/VR_method_comparison_2018/
plus ZIP. 42 CMS-style PNG/PDF pairs: 18 category/channel plots and 3
six-category overviews per calibration, two 21-page bundles, six chi2 tables.
Compare direct DY, direct pT-corrected DY, direct mass-corrected DY, nominal
alpha and pT-corrected alpha with signed Data-TTplusWW-VZ in VR. No ZH.

Direct VR numerical spectra were already in the channel-pT and fineAlpha-mass
archives, but not included in previous SR comparison plots. Validate and use
these rather than refit or scan trees. 36 category/channel/calibration input
matches (8 fields plus primitive hash) pass. All 180 chi2/rank values are
independently recomputed from residual covariance and populated-bin mask;
all 180 joint ratio variances are independently checked using observation,
prediction and residual covariance diagonals. Full numerical outputs packaged.
15 existing covariance/mass regression tests and script compilation pass.
Resolved-untag primary and combined CR-only overview visually inspected.

inclVR factors use C+VR: direct reweighting and pT-alpha VR are partly
in-sample, even though alpha itself is DY(VR)/DY(C) applied to Z(C).
heldoutVR factors use only C; prefer that transfer cross-check. Original
mass-fit adaptive boundaries used both calibration sets and therefore VR
statistics, so it is not a fully blind adaptive-grid validation. Bin-choice
uncertainty is excluded. Joint data-reuse covariance is retained in both.

Combined primary chi2/rank (DY / DYpt / DYmass / alpha / alphapt):
resolved btag:17.55/7,11.26/7,11.48/7,12.45/7,12.32/7;
resolved untag:441.70/10,180.34/10,169.80/10,167.24/10,159.74/10;
resolved VBF:15.18/9,16.98/9,15.85/9,17.54/9,16.37/9;
merged btag:9.93/3,7.13/3,4.85/3,4.52/3,6.17/3;
merged untag:15.84/4,5.66/4,6.46/4,6.59/4,6.47/4;
merged VBF:0.97/2,0.78/2,0.71/2,0.38/2,0.71/2.
CR-only resolved untag:441.70/10,174.38/10,167.23/10,167.24/10,156.34/10;
merged untag:15.84/4,5.22/4,6.46/4,6.59/4,6.27/4.
No universal winning method: modest improvements do not cure resolved untag
transfer nonclosure; merged pT-alpha can worsen btag. Statistical-only tests,
no VR normalization floated, original mass bins and full tails preserved.

## 2026-10-07: Reduced-tree region diagram, final blinded delivery

Standalone Plot/plot_region_definitions_reduced.py adapts the legacy
plot1DDataMCRegionALL without changing the legacy script/settings. Final
Plot/region_definitions_2018_blinded/ and ZIP contain 25 PNG/PDF pairs,
25-page bundle, histograms/covariances, bins/yields, provenance, tests and code.
Local mirror has the same folder name. Earlier region_definitions_2018 and
region_definitions_2018_final outputs are pre-instruction unblinded caches,
NOT delivery packages. Final SR70-105 Data is NaN in exported arrays/CSV,
hidden in both panels; Data-SR source-yield audit rows removed. No cache copy.
No per-bin vertical outlines or footer notes. Enlarged labels/ticks/legend;
region boundaries only, no105-135 removal hatch. Marked CR40-70, SR70-105,
CR105-135, VR135-165, CR165-180. CMS Work in progress,59.83fb^-1.

Current CR omits105-135 but includes VR; slice low/high current CR, bridge
105-135 from RedTrees_old broad CR, use current SR and VR separately. All48
sorted full (mass,mass2q,weight) common-CR multiset comparisons pass after
excluding bridge, m2l2q>=300. No event IDs: not an ID uniqueness audit. Source
slices are disjoint; nominal weights applied once; no reweighting, alpha or ZH.

Old get_err returns sqrt(hist variance): its total band is MC statistical,
not evaluated full systematics. Native revised band is sqrt(sum_process sumw2),
signed weights retained. Exact central68.27% Data Poisson intervals, MC-only
ratio band, Data-only point errors. Full systematic model unavailable; no
arbitrary percentages/residual inflation. Optional signed coherent process
normalization nuisances propagate full outer-product bin covariance.

All-category/channel aggregate sums entries from both topologies, not a
deduplicated event union: selection allows overlap, reduced trees lack IDs.
Per-process variance upper bound (sigma_resolved+sigma_merged)^2; process
bounds summed. Integrated bound uses integrated native variances, not the
sum of bin bounds. Archive exports variance envelope only for aggregate,
not fictitious off-diagonals. Aggregate Data sums central84.13% native Poisson
interval endpoints (Bonferroni >=68.27% coverage regardless of overlap).
Native plots retain ordinary Poisson intervals. Exact deduplication needs IDs.

Validation:10 unit tests,48 common-CR matches,25 histogram checks,180 audited
source records after blinding,25 PNG/PDF pairs and25 bundle pages. Aggregate
and resolved-all/combined layouts visually inspected: no SR points/ratios,
no bin-edge vertical lines, no footer, enlarged texts fit. Source code/tests
and comprehensive README retained remotely and in download. Original plotting
scripts/trees unchanged; additional validated intermediate outputs retained.

## 2026-10-07: Six 2018 CR nominal Data/MC category plots

Located matching screenshot routine Plot/plotMcData.py::PlotUnit.plot1DDataMC
(line107). Standalone adaptation Plot/plot_cr_datamc_reduced.py, tests and
CR_DATAMC_REDUCED.md retained; legacy routine unchanged. Final delivery is
Plot/cr_datamc_2018_CMS/ plus ZIP, mirrored locally, six individual PNG/PDF
pairs and all_CR_datamc.pdf six-page bundle. Earlier fullCR renders are layout
previews, not delivery. Six resolved/merged x btag/untag/vbftag plots combine
ee+mumu, nominal DY MC stacked with TTplusWW/VZ, no alpha/reweighting/refit/ZH.
Full CR40-70 plus105-180 includesVR. SR Data AND MC omitted from both panels,
tables and all histogram arrays. Reduced-tree range40-180, m2l2q>=300.

Uses audited region_definitions_2018_blinded cache, hash verified,48 bridge
common-CR matches retained. All24 combined-channel process histograms match
ee+mumu sumw/sumw2 exactly. No tree reads/production needed. MCstat sumw2
band and exact68.27% Data Poisson intervals; separate Data-only ratio errors.
Log upper panels, large text/no footer/bin outlines. Contiguous band segments
retain boundary bins besideSR. Ratio ranges retain sparse errors, no clipping.
ZeroData log markers denote upper endpoints, not invented central counts.
Six tests pass, six plot pairs plus bundle checked. Resolved-untag and merged-
VBF previews inspected; final two-column legend avoids overlap. Script/cache
are reusable for subsequent plots to avoid repeating discovery/preparation.

Latest user typography request: category annotation18->21pt, legend14->17pt.
Final delivery superseded by Plot/cr_datamc_2018_CMS_largeText/ and ZIP.
Style-only cached rerender; histogram/yield/bin values remain unchanged.

Presentation typography supersedes largeText: cr_datamc_2018_CMS_presentation_final/
plus ZIP. Category40pt, y44/ratio40/x46pt, ticks28pt: approximately1.8-1.9x.
Legend26pt (1.53x), CMS/luminosity30pt (1.36x) limited to avoid collisions.
Compact single-column legend, fewer x tick labels, added log headroom and
larger margins. Six individual PDF/PNG pairs and a vector3x2 six-panel page
plus six-page bundle. Style-only cached rendering; numerical tables unchanged.
Six-panel preview checked; final tight legend spacing avoids resolved-untag
Data overlap, and explicit ratio ticks retain the reference scale. CSV hashes
match largeText version. No additional physics/input validation scan needed.

Latest user request supersedes presentation_final: final delivery is
Plot/cr_datamc_2018_CMS_Preliminary/ plus ZIP. Category annotations40->30pt,
legend26pt and all axes/header text unchanged. CMS label changed to Preliminary
at user's explicit request. Individuals, six-page bundle and3x2 vector page
rerendered from same blinded cache; no numerical changes or new tree scans.

## 2026-10-07: Stored nominal 2018 alpha ratios, separate channels

New Plot/plot_stored_alpha_ratios.py and STORED_ALPHA_PLOTS.md. Directly read
only FitTemplate/AlphaRatio/AlphaHisto_SR_2018.root and AlphaHisto_VR_2018.root.
Keys resolved_isEE/isMuMu_TAG and merged_tag_isEE/isMuMu_TAG. No fixCR/oldRES,
2lep combined histograms, rebinning, fitting or recomputation. Central values
and explicit stored Sumw2 errors preserved. Files unchanged, digests audited.
Resolved37 bins, mergedEE20/mumu19 bins, stored range0-3500GeV. Display300-3500
search domain; all normal bins and exact edges exported, flow values audited.
Finite wide last-bin horizontal intervals retained, no invented overflow.

Delivery Plot/alpha_ratios_2018_CMS_Preliminary_search/ plus ZIP:24 individual
PDF/PNG pairs (SR/VR xee/mumu x6 categories), four vector3x2 category pages
and24-page bundle. Earlier full-range render is a preview, not delivery.
Latest CMS Preliminary style, category30/legend26/CMS30/y44/x46/tick28pt.
72 exported arrays validated exactly against source reads; stored statistical
errors are not a full systematic/cross-bin covariance model. No Data read,
no production ROOT outputs or jobs changed. Script uses unbin, reusable.

Legend-free delivery supersedes search package: Plot/alpha_ratios_2018_CMS_noLegend/
and ZIP. Removed redundant legends from all24 individuals and4 overview
pages; all other styling, nominal alpha values, edges and errors unchanged.
Caption notes now describe stored statistical y errors and original-bin x
spans. Same direct standard ROOT inputs, no rebinning/recomputation.

## 2026-10-07: Alpha-estimated SR/VR background stacks, combined channels

Delivery: Plot/alpha_background_2018_CMS_Preliminary/ and sibling ZIP.
Twelve individual PDF/PNG pairs, two vector six-category overview pages,
and a twelve-page PDF bundle. CMS Preliminary presentation styling retained.
Zjet is the sum of separately transferred ee and mumu source Data-TTplusWW-VZ;
target TTplusWW and VZ MC are added. No ZH subtraction or reweighting.
No Data_SR files opened; no target Data or ratio panels displayed.

SR alpha uses full CR including VR. VR alpha uses CR excluding VR.
Current low/high CR plus validated old-tree 105-135 bridge defines C;
SR source=C+VR, VR source=C. All48 common-CR observable/weight multiset
comparisons passed. Earlier nominal study outputs had coarser binning.
Recomputed DY ratios explicitly on original standard AlphaHisto_SR/VR_2018
channel-native edges, adjacent merging where source/target Neff<5 or signed
source nonpositive. Mass>=300 retained, >=3500 tail folded into last bin.
Original production ROOTs unchanged. Diagnostic ratios saved separately as
AlphaHisto_SR_fullCR_2018_diagnostic.root and
AlphaHisto_VR_CRnoVR_2018_diagnostic.root; not deployed to production.

Merged ee/mumu native edges differ: transfer each on its own native grid,
project to electron display edges and sum, retaining shared-alpha cross-bin
covariance. Exported NPZ contains native alpha edges and projected covariance.
Sparse merged VR btag has two displayed bins, vbftag one; not a shape claim.
Grey total statistical band includes source Data/MC subtraction, DY alpha MC,
and target TT/VZ MC statistics. Blue alpha-MC component is already contained
in grey total, not an additional band to add. No full systematic uncertainty.

Reusable Plot/plot_alpha_background.py, alpha_background_core.py,
test_alpha_background.py and ALPHA_BACKGROUND_PLOTS.md preserved.
Checks: seven regression tests,180 audited ROOT inputs,48 common-CR matches,
24 channel-transfer conservation checks,72 diagnostic ROOT roundtrip array
checks. Stored-alpha provenance comparisons exported rather than assumed.
Inputs/source histograms, predictions/covariances, yields/bins CSV, adaptive
binning JSON, hashes, metadata and validation report packaged with plots.

## 2026-10-07: VR observations and positivity-only plotting bins

Latest delivery supersedes the Neff>=5 preview:
Plot/alpha_background_2018_CMS_positive_VRdata/ and sibling ZIP.
Observed raw Data (ee+mumu) added to all six VR plots, with central 68.27%
Garwood Poisson intervals. SR Data remains unread and blinded. VR observations
do not enter VR source, prediction or bin choice. No ratio panel added.

By explicit user request --min-effective defaults to0: adjacent native bins
merge only when source DY, target DY or signed source Z is nonpositive.
Retain large statistical error bands rather than merging to Neff>=5.
Original conservative merged VBF VR single-bin result was caused by both
sparse target DY and source Z, compounded by greedy adjacent merging; it
was not a physical single-bin constraint. Positivity-only latest VR display
bins resolved untag/btag/VBF=31/30/31; merged=16/11/9. SR=31/29/31 and20/19/18.
Combined predictions still use independently transferred native-channel alpha
and full projected covariance. Full CR for SR; disjoint CRnoVR for VR.

Twelve individual PNG/PDF pairs, two six-category pages, twelve-page bundle,
diagnostic alpha ROOTs and numerical archives preserved. Observed counts
exported in histograms.npz and bins.csv; SR observed fields are empty.
Total VR observed counts resolved=44496/888/2783, merged=168/42/15.
Seven regression tests,48 CR matches,24 transfer conservation checks and72
diagnostic ROOT array roundtrips passed. Six-category VR PNG visually reviewed.
No production templates/cards modified. Scripts and README retained.

## 2026-10-07: VR Data/prediction ratio panels

Latest delivery Plot/alpha_background_2018_CMS_positive_VRratio/ and ZIP.
Cached-only renderer plot_alpha_background_ratio.py adds lower panels to six
individual VR PNG/PDF plots and the vector six-category VR overview.
Also all_VR_alpha_background_ratio.pdf (six pages); twelve-page SR/VR bundle
updated. Original SR individual/overview plots copied byte-for-byte unchanged,
hash-checked locally. SR stays blinded with no observed-data ratio.

Ratio=raw observed Data/(alpha-estimated Zjet+target TTplusWW+target VZ).
Point errors use central68.27% Garwood Data intervals divided by nominal
prediction. Prediction uncertainty is not added twice to points: grey relative
total statistical band centered on1 and contained blue alpha-MC component.
Zero observed counts appear at0; nonpositive denominators yield masked NaN.
Default axis0-2: nine central ratios above2 marked by upward triangles;
uncapped values/errors and masks saved in128-row VR_ratios.csv. --ratio-max
can enlarge display range. No per-bin vertical grid or extra footer added.

NPZ predictions/observations, native alpha factors, full covariance, binning,
yield and bin tables byte-for-byte unchanged versus positivity-only package.
No tree reads or production recomputation. Three focused tests passed for
separate numerator/denominator errors, zero-count intervals and denominator
masking. Updated six-category PNG visually reviewed; package downloaded.
