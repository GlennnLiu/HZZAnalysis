# NewScheme: reduced trees and analysis plots

Canonical script home on lxlogin:

`/cms/user/liugeliang/CMS/HZZAnalysis/NewScheme`

Reorganized on 2026-10-14. Only code and supporting configuration snapshots are
organized here. Existing ROOT files, histogram caches, figures, job logs and
original YAMLs stay in their previous locations. No production, submission,
refitting or plot regeneration is performed by the migration itself.

## Layout and entry points

| Directory | Main scripts and purpose |
|---|---|
| `production/background/` | `produce_reduced_stream.py` runs chunked DY/TTplusWW/VZ/Data reduction. `make_reduced_trees.py` implements cuts, variable mappings and weights. Frozen `cards/` and `sf/` are beside the scripts. |
| `production/signals/` | `prepare_signal.py` adds sample-specific DeepJet weights to cloned signal ntuples; `produce_signal_reduced.py` produces the signal reduced trees. Frozen cards/SFs and reference helpers are included. |
| `plots/z_kinematics/` | `plot_run2_z.py`: leptonic/hadronic Z pT and mass. |
| `plots/discriminants/` | `plot_run2_kd.py`: KD, jjVBF and jVBF discriminants. |
| `plots/pnet/` | `plot_run2_pnet.py`: selected ParticleNet scores, including DY SF-on/off comparison. |
| `plots/pnet_fullrange/` | `plot_fullrange_pnet.py`: pre-working-point ZvsQCD distribution from original ntuples, not selected reduced trees; no ParticleNet SF. |
| `plots/dy_pnet/` | `plot_dy_pnet.py`: DY distributions with/without ParticleNet, across eras. |
| `plots/sr_pnet/` | `check_sr_pnet.py`: direct SR data/MC ParticleNet comparison. This diagnostic uses SR observations; it is not the blinded presentation workflow. |
| `plots/alpha_pnet/` | `check_alpha_pnet.py`: alpha ratios with/without DY ParticleNet SF. |
| `plots/alpha_closure/` | `check_alpha_closure.py`: SR/VR data/background alpha closure. `replot_cached.py` redraws saved histograms. This diagnostic uses target observations. |
| `plots/signal_overlay/` | `prepare_histograms.py`: caches signal and SR mass spectra; `render_plots.py`: initial signal overlays; `revise_layouts.py`: latest in-frame layout and CR mass plots. |
| `studies/alpha_ratio/` | Older alpha/ZH subtraction, auxiliary CR and pT/mass reweighting studies. `redraw_without_mass.py` is the current four-method VR display, with optional enlarged individual text. |

Tests and validation helpers are retained beside each workflow. The mapping
from old to new paths and original hashes is in `SOURCE_MANIFEST.json`.
`MIGRATION_VALIDATION.json` records the relocation checks.

## Environment

```bash
source /scratchfs2/cms/liugeliang/miniconda3/bin/activate unbin
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
HZZ_ANALYSIS=/cms/user/liugeliang/CMS/HZZAnalysis
cd "$HZZ_ANALYSIS/NewScheme"
```

Run a primary script with `--help` to inspect its options. Do not execute
`submit_jobs.py` casually: it is a historical submission helper, protected by
the retained `submitted_jobs.json`, not a plotting command. Batch files are
historical reproduction references; their existing output destinations must
not be resubmitted. The migration does not submit jobs.

## Reduced-tree inputs/outputs and weights

Existing validated background/Data trees:
`$HZZ_ANALYSIS/FitTemplate/RedTrees_20261009_pnet/`.
Existing signal trees: `$HZZ_ANALYSIS/FitTemplate/RedTrees_signal_20261012/`.
Full DeepJet-augmented signals:
`$HZZ_ANALYSIS/FitTemplate/SignalWithDeepJet_20261012/`.

Background trees save `weight`, `weight_no_pnet`, `pnet_sf`,
`pnet_sf_zvsqcd` and `pnet_sf_zvslight`. `weight=weight_no_pnet*pnet_sf`.
Merged DY uses the frozen era/flavour/pT/working-point tables, not the legacy
constant 1.36. ZvsQCD is applied to merged categories; ZvsLight additionally
to merged non-VBF categories. Resolved ParticleNet factors are one. DY DeepJet
remains unapplied, matching that production. Negative MC weights are retained.

Signals are ggF/VBF at 500/1000/3000 GeV for all four eras, with separate
DeepJet augmentation and a 1 pb benchmark normalization per hypothesis.
DeepJet is applied only to resolved b-tag/untag, never the invalid sentinel
domain. See `production/signals/README.md` for the full weighting audit.

CR is the existing hadronic sideband outside strict 70–105 GeV, including
105–135 and VR. VR is strict 135–165 GeV. For held-out VR alpha, remove VR
from full CR; never add separate VR trees to a full-CR tree. All m2l2q>=300
belongs to the search spectrum, with no signal-free low-mass slice.

Example new background production (choose a genuinely fresh output):

```bash
python production/background/produce_reduced_stream.py \
  --era 2018 --cat merged --output-dir /scratchfs2/cms/liugeliang/NEW_BACKGROUND_REDUCTION
```

The runner now defaults to its local frozen `cards/` and `sf/`, independent of
the working directory. It refuses to overwrite final or partial trees.
Signal augmentation and reduction also require explicit output directories;
never reuse the existing production roots for a new run.

## Plot examples: newer ParticleNet-corrected trees

```bash
python plots/discriminants/plot_run2_kd.py \
  --input "$HZZ_ANALYSIS/FitTemplate/RedTrees_20261009_pnet" \
  --output /scratchfs2/cms/liugeliang/NEW_KD_PLOTS

python plots/pnet/plot_run2_pnet.py \
  --input "$HZZ_ANALYSIS/FitTemplate/RedTrees_20261009_pnet" \
  --output /scratchfs2/cms/liugeliang/NEW_PNET_PLOTS
```

These distribution workflows combine the four eras, tags and lepton channels,
display 138 fb^-1, blind SR and show CR data. Native invalid KD/score values
are excluded, never folded into physical bins. Selected ParticleNet plots
toggle DY only; other MC stays unchanged. Most plotters support `--histograms`
to redraw an existing cache without reading ROOT trees again.

For signal overlays, prepare or reuse the existing cache at
`$HZZ_ANALYSIS/Plot/signal_overlay_Run2_20261013/`. Example latest layouts:

```bash
python plots/signal_overlay/revise_layouts.py --analysis "$HZZ_ANALYSIS"
# Only CR mass data/MC, direct DY MC, no signals:
python plots/signal_overlay/revise_layouts.py --analysis "$HZZ_ANALYSIS" --part cr
```

These commands write into the existing revised `Plot/*_layout` folders and
`Plot/run2_CR_m2l2q_DataMC/`; they are not read-only. SR remains blinded,
signals appear only in SR, and CR has observed data and data/MC panels.
Mass bins are event-filled common bins across channels/tags, with the >=3500
tail retained and density-scaled to Events / 50 GeV. No bin interpolation is
used to combine channels. Legends remain inside the frame; measured text
clearance sets the range and CMS labels remain at the border.

## Older method-comparison study: do not confuse inputs

The saved 2018 pT-reweighting/method-comparison study predates the newer DY
ParticleNet-corrected production. Moving its scripts does not update its inputs
or turn it into a new ParticleNet-corrected study. Original cached results stay
under `$HZZ_ANALYSIS/FitTemplate/AlphaRatio_examination/plots/`.

Primary correction code: `studies/alpha_ratio/pt_reweight_channel.py` and
`pt_reweight_core.py`. Comparison plotters: `plot_sr_methods.py` and
`plot_vr_methods.py`. ZH statistical-band display: `plot_sr_zh_statband.py`.

Current four-method VR presentation, from saved results without refits:

```bash
python studies/alpha_ratio/redraw_without_mass.py \
  --input "$HZZ_ANALYSIS/FitTemplate/AlphaRatio_examination/plots/VR_method_comparison_2018" \
  --output /scratchfs2/cms/liugeliang/NEW_VR_DISPLAY --large-individuals
```

Four methods remain: direct DY, pT-reweighted DY, nominal alpha and
pT-reweighted alpha. `heldoutVR` fits factor values excluding VR;
`inclVR` fits include VR and are partly in-sample. In both, VR alpha source
excludes VR. Shared prediction/observation covariance is retained. The chi2
denominator is residual-covariance rank, not a fitted-parameter-adjusted ndf.

Without `--large-individuals`, the redraw creates all summary/individual
figures. With it, only individual files and bundles are updated; standalone
summary files stay unchanged. The output contains no raw ROOT files.

## Compatibility and preservation

Old lxlogin plot/study Python entry points forward to their NewScheme locations
through compatibility links. Their exact original sources are retained as
`.py.pre_NewScheme_20261014` backups; the installer refuses changed sources.
Frozen production job directories are deliberately untouched, preserving
historical code hashes, original YAMLs, submission records and logs.
New production work should use `NewScheme/production`, not those old snapshots.

No ROOT files, histogram caches or physics calibration values are moved or
deleted. Local old script copies remain reference snapshots; use the local
NewScheme mirror for maintained code. Compatibility details are recorded in
`COMPATIBILITY.json`. Source-only provenance is in `provenance/`.
