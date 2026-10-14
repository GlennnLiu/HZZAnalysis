# Layout review and CR mass plots — 2026-10-13

Inspected signal and previous non-signal plots. Fixed y-limit headroom factors of
100–10000 were squeezing histograms; legends also collided with category text
and high signal bins. An external-header preview was rejected by the user.
The final renderer keeps legends/category text inside the plot and CMS labels
at the border. Canvas-measured text boxes determine the required clearance,
using the actual yields and errors locally underneath each text box. The
legend is placed on the side requiring less y-range. Original figure sizes
and plotting borders are retained.

All revised plots are retained in separate `_layout` folders. Original plotting
scripts, caches, ROOT files and configuration YAMLs are unchanged. No new
ntuples, skims or batch jobs are produced. CR signal overlays are removed from
all revised versions; SR keeps the six independent 1 pb signal hypotheses.

New Run-2 CR mass plots: `../run2_CR_m2l2q_DataMC`. Existing
`RedTrees_20261009_pnet/*_CR.root` inputs are histogrammed once using nominal
weights, including DY ParticleNet. Selection: m2l2q>=300 and 40<=mass2q<180;
the stored SR exclusion is asserted. The full CR includes 105–135 and VR;
separate VR trees are not loaded. No SR data is read. Data counts use unit
weights and Poisson errors; MC uncertainty is sqrt(sumw2). Signed weights
remain intact. Mass bins match the previous signal-overlay display, with
common event-level edges and a retained >=3500 tail. Independent assertions
check source yield and squared-weight sums per topology/tag/channel/process.

24 individual CR mass plots cover both topologies, three tags plus inclusive
tags, and ee/mumu/combined channels. Three six-category overview pages and
one inclusive overview are provided, along with a 24-page PDF bundle. Every
figure has PNG and PDF versions and a data/MC panel. Reproduction:

```
python revise_layouts.py --analysis /cms/user/liugeliang/CMS/HZZAnalysis
python package_layouts.py
```

Use `--part signals`, `--part backgrounds`, or `--part cr` to update just one
group from cached inputs. Layout validation and output counts are retained in
JSON files alongside the code.

Completed: 13 output folders, 144 PNGs and 157 PDFs. Input yield/variance
assertions and measured in-frame clearance checks passed. Representative
signal/no-signal MELA, ParticleNet and CR mass PNGs were visually inspected.
