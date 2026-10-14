# Full-range Run-2 ParticleNet Z-vs-QCD

Read the original combined ANATree ntuples, not the working-point-filtered reduced trees. Four eras, ee+mumu and all three exclusive tags are combined. Only the `(particleNetZvsQCD>0.9)` requirement is removed from each original merged tag/channel cut. All other era-specific cuts are retained, plus m2l2q>=300 and 40<=massmerged<180. Each event is counted once; tag/channel overlaps are checked. CR includes the 105–135 bridge and VR; SR is blinded.

No ParticleNet SF is applied to any MC in this full-range baseline. Use EventWeight*lumi*1000*xsec/sumWeights from exactly one combined file per sample. Existing lepton/event corrections in EventWeight are retained. No fail-region SF is guessed or extrapolated. Original luminosity normalization is retained and the displayed label is 138 fb^-1. Stack colors and text conventions match preceding plots. The band is MC statistical uncertainty sqrt(sumw2), not systematics.

The original Data ntuple contains all regions. Only CR cut expressions are evaluated for Data, and only CR events enter histogramming/audits. No SR Data histogram, selected yield, distribution or ratio is produced. SR MC is shown without observed data.

Twenty bins covering [0,1], including exact endpoints; invalid scores are excluded and audited. Validation checks full and tagged weighted yields and variances. The score>0.9 subset is compared bin-by-bin with prior cached no-PNet DY distributions for SR/CR and with prior CR Data exactly. This is a selection/normalization regression test, not an SF test.

Existing reduced trees and YAMLs are read-only and unchanged. No new reduced-tree production is performed. Exact input filenames and YAML text/hash snapshots are stored in metadata.json; histogram cache and per-sample yield audit are retained.

Outputs: separate `merged_{SR,CR}_particleNetZvsQCD_fullrange` PNG/PDF pairs; a `SR_CR_fullrange_particleNetZvsQCD` overview pair; two-page `all_Run2_fullrange_PNet.pdf`.

```bash
/scratchfs2/cms/liugeliang/miniconda3/envs/unbin/bin/python plot_fullrange_pnet.py --input /cms/user/liugeliang/CMS/HZZAnalysis --output .
```

Add `--histograms histograms.npz` to rerender from cache without ROOT reads.
