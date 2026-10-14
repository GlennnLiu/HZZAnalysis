# Run-2 leptonic/hadronic Z kinematics

Four eras combined: 2016preAPV, 2016postAPV, 2017 and 2018. Combine ee+mumu and
the three mutually exclusive main tagging categories (untag, btag, vbftag),
keeping resolved and merged topologies separate. Inputs are the validated
`FitTemplate/RedTrees_20261009_pnet` trees; no new production is needed.

All stacks use direct MC: DY+jets (green), TT/WW/tZq (purple), VZ=WZ+ZZ (orange).
The full nominal `weight` includes the approved era-specific ParticleNet SF for
DY merged events and the existing applicable corrections for other MC. No alpha
method, DY pT/mass reweighting or ZH subtraction is applied. Signed MC weights
are retained; the band is sqrt(sumw2), statistical only.

Selection retains original era/category cuts and requires m2l2q>=300 GeV and
40<=m2q<180 GeV. SR is strict 70<m2q<105, and CR is m2q<70 or m2q>105,
including 105–135 and VR135–165. The VR trees are NOT added to CR, which would
double-count them. Existing older-era cut differences are retained, not silently
harmonized. Normalization uses each era's original luminosity; their sum is
137.64 fb^-1, displayed as **138 fb^-1** as requested. No additional rescaling.

SR is blinded: SR Data trees are never opened, there are no SR Data histogram
keys, and no SR Data counts or yields enter the plots/reports. SR figures show
MC only, without a Data/MC ratio. CR figures show Data and the full MC stack,
with a ratio panel, using Garwood Poisson errors for Data and a separate MC
statistical band around one. Ratios zoom to 0–2.5, with upward triangles marking
off-scale central ratios. Full values remain in the histogram archive.

## Plots

Four variables: `pt2l`, `mass2l`, `pt2q`, `mass2q`. The hadronic quantities are
dijet pT/mass in resolved, AK8-candidate pT/mass in merged. Sixteen separate
PNG/PDF pairs: `{resolved,merged}_{SR,CR}_{variable}`. Four overview pairs:
`{resolved,merged}_{SR,CR}_four_variables`. `all_Run2_Z_kinematics.pdf` contains
the sixteen individual figures. CMS Preliminary styling and previous stack colors.

pT bins are 50 GeV through 1000, then 1000–1200, 1200–1500, 1500–2000 and
>=2000. The open final bin is displayed as 2000–2500, labelled >=2000; every
higher-pT event is retained. Plot density is Events/50 GeV. Leptonic mass bins
are 2 GeV over 40–140. Hadronic mass bins are 5 GeV over SR70–105 or CR40–180;
the excluded 70–105 interval is left empty in the CR plot for both Data and MC.
Variable-width bins are density-scaled; archives retain unscaled yields/variances.

`histograms.npz` stores edges and (sumw,sumw2) arrays. `input_yields.csv` stores
per-input yields and CR bridge/VR membership. `metadata.json` records the exact
336 input files and hashes of eight production validation reports, including
the explicit SR_Data_read=false flag. `plot_run2_z.py` retains reproducible code.

## Reproduction

```bash
cd /cms/user/liugeliang/CMS/HZZAnalysis/Plot/run2_z_kinematics_PNet
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python plot_run2_z.py \
 --input /cms/user/liugeliang/CMS/HZZAnalysis/FitTemplate/RedTrees_20261009_pnet \
 --output /path/to/fresh/plots
```

For style-only updates, add `--histograms histograms.npz`; this reuses the saved
bin sums without ROOT reads. No old plotting packages are changed.
