# Run-2 signal overlays — 2026-10-13

Requested signal additions are saved separately; previous plots, reduced trees and configuration files are unchanged.

- Six hypotheses: ggF and VBF at 500, 1000 and 3000 GeV, each normalized to a 1 pb benchmark. Dotted ggF and solid VBF share the mass colors red, dark red and pink.
- Four eras are combined using their original luminosity normalizations; the plot luminosity is displayed as 138 fb^-1.
- MELA and selected ParticleNet backgrounds/data are reused from validated histogram caches. Signals use nominal corrected reduced-tree weights. Original CR data/background ratios and statistical bands are unchanged; signals are not included in the ratio denominator.
- Selected ParticleNet plots retain the DY-SF-off total comparison. Full-range ZvsQCD plots remove only its >0.9 selection and use no ParticleNet SF, matching the previous full-range convention. Both logarithmic and linear variants are provided.
- New SR mass plots use direct DY MC, as explicitly confirmed by the user, with TTplusWW and VZ. SR data is never loaded or drawn. Individual channels/tags, combined channels, inclusive tags and overview pages are supplied.
- Mass plots use common event-filled edges: 50 GeV bins from 300 to 1500, 100 GeV bins to 3500, and an open overflow >=3500 displayed in 3500–4000. Histograms are density-scaled to Events / 50 GeV. This grid resolves the 3000 GeV hypotheses; it does not interpolate differently binned electron/muon histograms.

`prepare_histograms.py` builds caches, input yields and source provenance. It asserts full-range/tagged signal consistency, exact mass yield/variance accounting and SR blinding. `render_plots.py` reuses the caches for fast style changes. The three retained legacy plotting modules preserve previous backgrounds, colors, validity masks and CR ratios.

Production inputs: `FitTemplate/RedTrees_signal_20261012`, `FitTemplate/SignalWithDeepJet_20261012`, and `FitTemplate/RedTrees_20261009_pnet`. All six signal hypotheses are available for all four eras. No new batch production is required.

Histogram production passed all numerical assertions. Rendering outputs and plot counts are recorded in `outputs.json`; `metadata.json` and `input_yields.csv` retain the source audit.
