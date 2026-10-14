# Full-range ParticleNet plotting

Read original combined ANATree ntuples with the frozen reduced-production cut YAMLs. Followed the streaming producer's explicit 11-sample inventory without editing any YAML. Removed only the Z-vs-QCD >0.9 term. Combined four eras, both channels and all exclusive tag cuts; checked no event enters overlapping tags/channels.

No ParticleNet SF is applied anywhere in this baseline, because the current Z-vs-QCD calibration helper has no fail-region coverage. Kept EventWeight corrections and original cross-section/luminosity normalization. Accumulated weights and sumw2 in double precision. Used two bounded workers and retained per-era histogram caches.

Full-range yield/variance checks passed. Raw DY events with score>0.9 reproduce the previous cached no-PNet DY bin contents and variances in SR and CR. CR Data score>0.9 reproduces the previous observed histogram exactly. SR Data selections/histograms/yields are absent. Original Data ntuples contain all regions, but only CR cut expressions are evaluated for Data.

Existing trees and YAMLs were not overwritten. Saved only this isolated plotting package, histogram caches, source/config provenance and audits. No new ROOT/reduced-tree production was necessary.
