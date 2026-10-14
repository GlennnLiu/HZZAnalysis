# Mass-specific signal trees, 2026-10-12

Scope: ggF (`ggh`) and VBF, nominal masses 500/1000/3000 GeV, all four eras. Input files are Zhuhao's individual `*_M<mass>_*.root` skims, not narrow selections on `GEN_H1_mass`. Generated mass fluctuates; it is retained as `gen_mass`, separate from `nominal_mass`. Some combined-file components share SumWeight, so SumWeight alone is not used to split combined samples.

Original ntuples, shared helpers, background/data trees and original YAMLs are untouched. Code, cut YAMLs and ParticleNet YAMLs are frozen here. Output roots:

- `FitTemplate/SignalWithDeepJet_20261012/<era>/<ggh|vbf><mass>.root`: full cloned signal skims plus four DeepJet branches.
- `FitTemplate/RedTrees_signal_20261012/`: separate channel/topology/tag/SR-CR-VR reduced trees.

## DeepJet

Use each era's official medium-WP correctionlib `btagging.json.gz`: `deepJet_comb` for b/c and `deepJet_incl` for light jets. WPs are 0.2598/0.2489/0.3040/0.2783 for preAPV/postAPV/2017/2018. SF evaluation preserves the input jet/event order; it does not concatenate heavy/light selections as the legacy helper does.

Nominal resolved/non-VBF domain: Met<=150, found leptonic and dijet Z candidates, ee or mumu, KD_jjVBF<0.5 and KD_jVBF<0.7. For each mass-specific signal sample, measure unweighted passing/total jet counts by flavour, pT and |eta| before the b-tag cut, matching the previous ZH augmentation approach. Save the counts and efficiency grids. Do not borrow top/VZ/ggF efficiencies for another signal sample, clip weights or silently fill empty efficiency bins.

Per-jet weight is SF for a passing jet; for a failing jet it is `(1-SF*epsilon)/(1-epsilon)`, where epsilon is the passing probability. Multiply the two jet weights. Save `w_deepjet`, `w_deepjet_up`, `w_deepjet_dn` and `w_deepjet_valid`. Outside the domain, raw weight branches have -999 sentinel and validity is false. Apply DeepJet only to resolved b-tag/untag reductions. The applied factor is explicitly 1 for resolved VBF-tag and all merged events. Up/down branches vary calibration SFs, not the finite-MC efficiency uncertainty; efficiency counts are retained for that separate study.

Regression against the corrected ggF/VBF combined skims tests nominal/up/down two-passing-jet weights for every era, independent of efficiency assumptions. Untagged weights may differ from historical combined weights because these use mass-specific, correctly aligned efficiencies. The shared legacy helper and stored combined weights are not modified.

The augmented output is a ROOT fast clone, including unrelated branches and metadata objects; only the four new DeepJet branches are appended. Entry count/order, IDs, normalization and key physics branches are validated against the original, and the appended weights are fully read back.

## Reduced trees

Reuse the prior reduced-tree variable aliases and frozen SR/CR/VR cuts, with explicit channel guards. CR includes 105–135 and VR; VR is 135–165, SR 70–105 GeV. Do not add VR trees to full CR. No alpha correction is applied to signal.

Merged signal: ParticleNet ZvsQCD for all tags; ZvsLight additionally for b-tag/untag, using the same frozen era/flavour/pT/WP helper as background production. Resolved ParticleNet factor is 1.

`weight_base = EventWeight * lumi_fb * 1000 * xsec_pb / SumWeight`, using each event's component SumWeight, not an aggregate TH1. Mass-specific signal card entries use **1 pb benchmark normalization**, not a claimed physical signal cross section. `weight_no_pnet=weight_base*deepjet_sf`; `weight=weight_no_pnet*pnet_sf`; `weight_no_deepjet=weight_base*pnet_sf`. Save all weight factors and raw DeepJet variation branches, `SumWeight`, nominal/generated mass and source entry index.

All original requested physics branches are retained: m2l2q, leptonic/hadronic Z pT and mass, kd, kd_jjvbf, kd_jvbf and both ParticleNet scores. Source KD sentinels are retained, not folded into physical values.

Expected production: 24 augmented skims, 864 reduced trees (108 per era/topology job), eight validation reports and completion markers. Each output is fully reread, checked for region leakage, finite values, weight-factor consistency, source-entry order and weighted-yield agreement. Partial files are only renamed after validation; existing outputs are never overwritten.

## Run

Use the `unbin` environment and one numerical-library thread. First run `prepare_signal.py --reference-only`, then `--pilot`; test both reduced topologies with `--sample ggh500`. Full augmentation uses two bounded workers. Only after all 24 augmented inputs are validated, submit one `run_<era>_<topology>.sh` per era/topology with hep_sub, 6000 MB. Preserve job IDs, logs and independent final-validation reports.
