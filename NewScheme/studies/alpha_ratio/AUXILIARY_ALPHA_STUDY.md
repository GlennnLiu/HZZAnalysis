# Independent-selection alpha closure: 2018 untagged

## Definition

The user's intended double ratio is measured in an orthogonal event selection
A, not in a low m(2l2q) control interval. For each total-mass bin x, define
Z = Data - TTplusWW - VZ (TTplusWW includes WW and tZq), and

```
alpha_A(x) = DY_A,target(x) / DY_A,sidebands(x)
kappa_A(x) = Z_A,target(x) / [alpha_A(x) Z_A,sidebands(x)]
           = [Z_A,target / Z_A,sidebands] / [DY_A,target / DY_A,sidebands]
candidate nominal prediction = kappa_A(x) alpha_N(x) Z_N,sidebands(x)
```

The essential additional assumption is kappa_A(x) = kappa_N(x). Observing
nonclosure in A does not establish this. Nominal target data are used here
only to diagnose that assumption, not to fit the auxiliary correction.
Nor does an orthogonal selection automatically establish negligible signal
leakage: a calibrated auxiliary measurement must assess that separately.

The nominal selection N is MET <=150; A is MET >150. All other current
same-flavour dilepton, non-VBF and untagged requirements are kept. Auxiliary
slices are 150 < MET <=250 and MET >250. Additional nominal MET slices are
diagnostic only, **not signal-free calibration regions**.

All m(2l2q) >=300 GeV is retained, including an explicit >=3500 bin.
The primary independent hadronic sidebands are 40-70,105-135,165-180 GeV;
target SR is 70-105 and independent target VR is 135-165. CSV files also
include sidebands without 105-135, and SR tests including the VR in the
source; the overlapping source/VR closure is deliberately not computed.
Histogram interval boundaries follow the left-closed convention.

## Candidate-flag audit

The original skim scan checked all eleven enabled samples. None of the
`foundZ1LCandidate == False` events had `isEE` or `isMuMu` true. The source
`ANATree/src/TreeLoop.cc` permits opposite-sign different-flavour pairs,
then sets `foundZ1LCandidate` only for same-flavour opposite-sign pairs;
the different-flavour pairs are marked `foundTTCRCandidate`.

Thus flipping this flag with the same channel requirements gives an empty
sample, not a DY-enriched alternate Z selection. Its natural role is a top
control measurement. Using it for that purpose also requires checking the
validity of category discriminants and weights for the different-flavour
events; it must not be treated as a drop-in copy of the nominal DY channel.

## DeepJet issue and status

The input high-MET resolved non-DY entries have `w_deepjet=-999` sentinels.
`btageff.py::SFs.loop` explicitly computes the branch only for MET <=150.
Using the stored branch after inverting MET therefore produces invalid
background subtraction. No input ROOT trees were modified.

The suggested `lib/deepjet.py` was inspected. Its `get_deepjet_sf` evaluates
heavy and light flavours separately and then concatenates the results;
assigning those concatenated values back to the original mixed-flavour
array does not preserve event order. Its untag helpers also use
`eff_fail = 1 - histogram_content`, whereas the stored-weight producer
uses the passing efficiency directly in (1-SF*epsilon)/(1-epsilon).
The efficiency producer defines its histograms as pass/total. These
differences must be resolved before directly using that helper for a
calibrated correction. The helper itself was not changed.

A diagnostic vectorized extension reproduces the actual stored-weight
producer, including its legacy both-fail second-jet charm-to-b lookup.
It was validated against 152169 nominal TT entries (maximum relative
difference 5.94e-8) and 331162 nominal WZ entries (5.96e-8). However it stops
on some selected ZZ events whose fail-jet efficiency lookup is exactly 1.
The nominal-selection efficiency calibration is insufficient for a safe
automatic extension into this auxiliary selection. No arbitrary clipping
or substitution is permitted in the strict mode.

**The saved resolved auxiliary results are therefore an explicitly
uncalibrated feasibility study:** high-MET resolved non-DY DeepJet SFs are
set to 1; nominal entries keep the original valid weights. Plot titles
and `source_audit.json` state this. Merged weights follow the existing
ParticleNet prescription and do not use DeepJet. All errors below are
statistical only; subtraction-scale variations are stress tests, not
official systematic uncertainties.

## Preliminary feasibility result

Combined ee+mumu, resolved untagged, all total masses >=300:

| Selection | Target | Integrated kappa | Statistical error |
|---|---|---:|---:|
| nominal MET <=150 | SR | 0.9734 | 0.0064 |
| auxiliary MET >150, unit DeepJet SF | SR | 1.1501 | 0.3018 |
| nominal MET <=150 | VR | 0.9204 | 0.0064 |
| auxiliary MET >150, unit DeepJet SF | VR | 0.5901 | 0.2003 |

These integrated values use the ratio of integrated yields, NOT the sum
of mass-bin predictions. The mass-bin calculations and separate channels
are all in `alpha_closure.csv` and the plots. Neither type is an application
of the stored production alpha histogram: these are newly computed
coarse-bin nominal DY ratios, for consistent selection comparisons.

The auxiliary SR has Data=150, TTplusWW=82.85 and VZ=1.27, leaving
Z=65.88 +/-12.38; its sidebands have Data=365, TTplusWW=223.46 and VZ=2.73,
leaving Z=138.81 +/-19.32. The corresponding DY yields are only 25.28 in
SR and 61.26 in the sidebands. The integrated alpha prediction is
57.28 +/-10.49. Thus top subtraction dominates the control-region data
composition, and the correction uncertainty is about 26% even before
calibration and background-normalization uncertainties. In individual
mass bins the precision is substantially worse, with negative signed
targets and invalid signed source denominators retained/flagged.

Changing TTplusWW normalization to 0.8/1.2 moves the auxiliary integrated
SR kappa to 1.089/1.270. These are deliberately illustrative stress tests.
The auxiliary result neither establishes matching nonclosure nor supports
a few-percent nominal correction. A nominal-versus-auxiliary compatibility
plot with large error bars is not evidence for transferability.

## Recommendation

The independent-selection double-ratio construction is sensible, but
this MET-inverted sample cannot currently determine a useful mass-dependent
few-percent correction. First resolve the auxiliary SF calibration and
constrain top subtraction with an appropriate different-flavour control
measurement. Then test transferability across MET slices and channels,
including signal leakage and any additional backgrounds relevant to the
altered selection. A simultaneous likelihood would preserve the sizable
auxiliary statistical and nuisance correlations better than treating its
central-value double ratio as an exact multiplier.

If that still leaves insufficient information, retain the nominal alpha
prediction with a validated closure/transfer uncertainty; do not inject
these noisy auxiliary central values into the production shape. An
alternative auxiliary selection must both suppress signal and retain
enough DY with similar hadronic-response correlations. This dataset audit
does not yet establish such an alternative; no total-mass interval above
300 GeV is being labelled a signal-free CR.

## Artifacts and reproduction

All outputs are under `plots/auxiliary_MET_untag_2018/`, both on lxlogin
and in the local mirror. PDF and PNG are provided for each channel/category:

- `resolved_combined_SR_MET_transfer.pdf`: nominal and auxiliary closure
  factors, then the nominal/auxiliary comparison.
- `resolved_combined_SR_auxiliary_closure.pdf`: signed target versus alpha
  prediction in nominal/auxiliary selections, and the provisional
  auxiliary-corrected nominal prediction (diagnostic only).
- `yields.csv`, `alpha_closure.csv`, `subtraction_stress.csv`.
- `source_histograms.npz`, `source_audit.json`: compact sumw/sumw2 cells,
  source paths, entry counts, normalization, configuration and weight status.

Five deterministic unit tests passed for signed subtraction, variance
propagation, invalid denominators and the MET boundary. Input loading uses
combined sample files only, chunked reads and two decompression threads.

```
source /scratchfs2/cms/liugeliang/miniconda3/bin/activate
conda activate unbin
cd /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/AlphaRatio_examination
python -m unittest test_auxiliary_alpha_closure -v
# From scratch, explicitly provisional:
python auxiliary_alpha_closure.py --auxiliary-unit-deepjet \
  --output plots/auxiliary_MET_untag_2018
# Replot saved cells; explicit reuse, no freshness check:
python auxiliary_alpha_closure.py --reuse-cache \
  --output plots/auxiliary_MET_untag_2018
```

Without `--auxiliary-unit-deepjet`, strict in-memory DeepJet recomputation
is attempted and guarded by agreement with stored nominal weights; with
the current efficiency files it fails as described above. No production
templates or skims were regenerated, and no batch jobs were submitted.
