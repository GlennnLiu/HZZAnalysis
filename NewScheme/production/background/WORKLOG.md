# Reduced-tree production worklog

## 2026-10-09 — all-era ParticleNet production

User approved applying the existing era/flavour/pT/working-point ParticleNet
calibration to DY and requested full and ParticleNet-excluded weights, plus
reduced trees for every era. No constant DY 1.36 factor is used in this run.

Implemented an isolated version of `make_reduced_trees.py` and a bounded-memory
`produce_reduced_stream.py`. Production explicitly loads all five DY pT bins,
TTTo2L2Nu, WWTo2L2Nu, tZq, WZTo2Q2L, ZZTo2Q2L and Data, even where an older
card's Samples_lists includes only DY. Combined skims only; signed MC weights
and the original normalization are preserved.

Added 11 physics fields and five weight fields, including `weight_no_pnet` and
separate ParticleNet factors. Both weights include any other existing corrections.
Added strict calibration coverage and nonfinite-weight checks. Explicit channel
guards prevent older card omissions from contaminating ee with mumu. Other
era-specific cuts and the existing DeepJet implementation are not changed.

Corrected older-era VR selections to 135–165 GeV in three cards; original cards
have `.before_VRfix_20261009` remote backups. CR retains 105–135 and VR. No new
total-mass signal-free control region, range cut, alpha files or plots are produced.

Validation: five unit tests passed. A 2018 merged pilot with 3,000 entries/sample
produced 72 trees, with 613 selected rows across its output slots; full reread and
schema/region/weight checks passed. Production code, four cards and four SF tables
are frozen in `FitTemplate/jobs/reduced_20261009_pnet`.

All eight era/topology preflights passed: 88 combined-input header checks, all
requested branches present, valid normalizations, both SF tables available in
each era, and CR/VR definitions consistent with the requested conventions.

Eight production jobs submitted, cluster IDs 81981958–81981965. Job-to-era mapping
is in README.md. At the latest scheduler check all eight were idle, not failed or
completed. New output directory: `FitTemplate/RedTrees_20261009_pnet`; existing
reduced trees remain untouched. Expected 576 ROOT trees. Completion requires
eight validation reports and `.done` markers, followed by `validate_production.py`.

Follow-up status: all eight jobs finished, each with a passing full-output reread.
The read-only aggregate validator passed for all 576 final ROOT trees; no partial
files remained. Production is complete, superseding the earlier queued status.

The canonical legacy producer is not overwritten: run the retained streaming
entry point for this schema/production. The helper snapshot is imported by it;
legacy alpha-output functionality is not used. Complete-production status must
not be inferred from the passing pilot.
