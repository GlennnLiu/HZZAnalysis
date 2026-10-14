# 2018 ParticleNet alpha comparison — 2026-10-09

Recomputed alpha with/without DY ParticleNet on identical finalized reduced trees,
not by comparing old stored alpha values with new trees. SR uses full CR; VR
uses CR with strict hadronic 135–165 excluded event-by-event. 105–135 retained.
Hadronic acceptance 40–180, m2l2q >=300, channels independent. Only DY read;
neither observed Data nor non-DY backgrounds enter alpha.

Original alpha ROOT files supply their channel/category-dependent edges. Resolved
has a 300 edge, but merged does not: the original first search-overlapping bin
is intersected with the 300-GeV lower cut. All subsequent edges are preserved,
without merging. Graphs stop at original finite upper edges; integrated summaries
and saved total primitives retain all masses >=300, including beyond graph range.

Three covariance tests passed. Paired target/source weight covariances are included
in both alpha comparison ratio and difference; target/source regions are disjoint.
Resolved correction is exactly one and alpha is unchanged. Statistical errors only;
no SF-calibration/modeling systematics. Original alpha files remain unchanged.

Integrated merged alpha relative shifts (percent; ee / mumu):

| Tag | SR / full CR | VR / CR excluding VR |
| --- | ---: | ---: |
| untag | -0.3072 / +0.1707 | -1.5669 / -1.6198 |
| btag | +0.2074 / -0.1365 | +1.7410 / +3.1825 |
| VBF | -0.6282 / -0.0952 | +0.7913 / -0.4184 |

Integrated combined-DY diagnostic SR shifts: -0.0364%, +0.0234%, -0.3598%;
VR: -1.5940%, +2.6384%, +0.1249%. These combined ratios are not used to
transfer a combined CR prediction; channel alphas remain separate.

In individual original bins with DY Neff>=5 in both regions and both SF choices,
max absolute shifts in SR are 6.89/4.35% untag, 6.94/9.50% btag, and
4.77/6.07% VBF (ee/mumu). Corresponding VR maxima: 7.56/9.66%, 9.95/8.62%,
6.04/3.50%. All positive-yield bins are plotted, including less precise sparse
ones. Integrated cancellation therefore does not establish binwise invariance.

Saved 24 individual PNG/PDF comparisons, four overview pairs, two twelve-page
PDF bundles, per-bin and integrated tables, input alpha-file hashes, primitive
histogram archive, code, tests and read-only package validation.
