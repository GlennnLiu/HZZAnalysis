# 2018 untagged alpha closure study

Script: `closure_untag.py`. Output folders: `plots/closure_untag_2018_current`
and `plots/closure_untag_2018_broadCR`. Every figure is saved as PNG and PDF.
`yields.csv` contains integrated results; `bins.csv` contains bin values and
variances; `input_audit.txt` records data-tree mass-window membership.

## Reproduction on lxlogin

```bash
source /scratchfs2/cms/liugeliang/miniconda3/bin/activate
conda activate unbin
cd /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate
python AlphaRatio_examination/closure_untag.py --base RedTrees \
  --output AlphaRatio_examination/plots/closure_untag_2018_current
python AlphaRatio_examination/closure_untag.py --base RedTrees \
  --cr-base RedTrees_old \
  --output AlphaRatio_examination/plots/closure_untag_2018_broadCR
```

## Definitions and checks

- Signed observation is Data minus TTplusWW minus VZ, retaining signed weights.
- Target SR: 70 < mass2q < 105 GeV; target VR: 135 < mass2q < 165 GeV.
- CR always excludes the target SR and VR. Broad CR retains 105--135 GeV;
  narrow CR removes it. Input mass2q acceptance is 40--180 GeV.
- Current `RedTrees` CR already removes 105--135 GeV. Therefore the two cuts
  give identical results there. The broad study uses `RedTrees_old` CR and
  verified current `RedTrees` SR/VR. Narrow-cut predictions match current
  predictions numerically for all tested category/channel/target combinations.
- Alpha is nominal DY(target)/DY(CR), recomputed in the same bins as the
  comparison. This is a controlled closure study, not an exact reproduction
  of the stored fine-bin production alpha. Ratio-after-rebinning and applying
  fine-bin alpha-before-rebinning need not be equal.
- Mass edges: 300,400,500,600,750,900,1100,1400,1800,2300,3500 GeV.
  Reported yields exclude flow bins. Upper mass edge is effectively exclusive
  apart from NumPy's final-edge histogram convention.
- Source trees are trusted to have correctly normalized MC and category SFs.
  This study does not independently validate the production weighting.
- Statistical variance of signed data is sum(w^2) for all three processes.
  Ratio variance uses independent DY numerator/denominator; prediction variance
  includes signed CR statistics and finite DY MC statistics. Statistical
  uncertainties exclude luminosity, SFs, TT/VZ normalization/shape, DY modeling,
  and correlations from systematic nuisance parameters.
- Chi-square values are approximate Gaussian diagnostics, especially unreliable
  in sparse high-mass or merged VR bins. They are not a final likelihood test.

## Findings

| Category | Channel | Target | Signed observed | Broad alpha prediction | Observed/predicted broad | Observed/predicted narrow |
|---|---|---|---:|---:|---:|---:|
| resolved | ee | SR | 21903.8 | 22384.7 | 0.9785 | 0.9241 |
| resolved | mumu | SR | 29032.0 | 29844.5 | 0.9728 | 0.9240 |
| resolved | ee | VR | 18182.3 | 19771.8 | 0.9196 | 0.8768 |
| resolved | mumu | VR | 24100.6 | 26336.9 | 0.9151 | 0.8739 |
| merged | ee | SR | 628.4 | 626.1 | 1.0036 | 1.0115 |
| merged | mumu | SR | 846.5 | 867.1 | 0.9762 | 0.9866 |
| merged | ee | VR | 68.9 | 52.1 | 1.3232 | 1.3112 |
| merged | mumu | VR | 86.6 | 75.0 | 1.1550 | 1.1771 |

Resolved broad-CR statistical chi-square / 10 bins: ee SR 39.9, mumu SR
44.7, ee VR 68.7, mumu VR 103.1. Direct DY values are respectively 252.7,
75.3, 371.4, 127.5. Improvement is clear for the broad CR, but not universal
for the narrow CR (muon cases worsen). Comparisons also differ in statistical
uncertainty because alpha adds CR and denominator MC statistics.

The two resolved lepton channels agree closely on integrated residuals. SR
closure changes with total mass; a simple normalization cannot fix the shape.
Merged SR yield closure is good; merged VR statistics are too limited to
establish a comparably precise modeling statement.

## Interpretation and improvements

If the true DY yield in region R and total-mass bin i is
`k_R(i) * DY_MC_R(i)`, then observed/predicted closure is
`k_target(i) / k_CR(i)` under exact other-background subtraction. Alpha cancels
a common correction but not a region-dependent discrepancy. The resolved
pattern is consistent with hadronic-mass / total-mass correlation mismodeling;
it does not identify the generator as the unique cause. Jet response,
selection efficiencies, flavor composition, and subtraction errors are also
possible. The missing SM ZH contribution in untagged categories is far too
small in the earlier audit to naturally explain residuals of this magnitude.

Priorities:
1. Retain 105--135 in the untagged CR: it improves precision and closure here.
2. Audit input normalization/SFs and TT/VZ variations before attributing the
   remaining residual to DY alone.
3. Measure signed data / DY across hadronic-mass slices and versus total mass,
   pT(Z), and jet kinematics. Distinguish normalization from shape changes.
4. Use a low-dimensional correction in correlated variables, or jointly fit
   smooth DY corrections across hadronic-mass sidebands, with appropriate
   uncertainty for interpolation into SR. Validate on held-out slices/channel
   or eras; use higher-mass VR primarily as extrapolation/closure validation.
5. A total-mass-only reweighting cancels from the binwise alpha ratio; improving
   transfer requires correcting correlated variables such as hadronic mass
   and pT. An empirical VR correction cannot automatically be copied into SR:
   their residuals here differ (roughly 8% vs 2--3%).
6. Avoid treating noisy fine-bin ratios as exact corrections. Smooth or fit
   with statistical regularization; retain transfer-model uncertainty.
7. SR data can be used for this diagnostic, but a correction fitted to SR
   can absorb a signal. Signal-injection tests or a simultaneous likelihood
   with explicit signal are needed before deploying an SR-derived correction.
