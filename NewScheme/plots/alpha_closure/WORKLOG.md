# Alpha-method compatibility with DY ParticleNet on/off — 2026-10-09

User explicitly requested observed SR and VR comparisons. Implemented
`B=alpha*(Data-TTplusWW-VZ)_source+(TTplusWW+VZ)_target`, varying only DY PNet
through alpha. Reused audited paired-alpha archive on original separate ee/mumu
grids. Source for SR is full CR; source for VR excludes VR event-by-event.
105–135 retained; 40–180 acceptance, m2l2q>=300. No ZH subtraction, new fits,
tree production or alpha deployment. Non-DY ParticleNet weights stay nominal.

Prediction variance includes source Data statistics, source MC subtraction,
DY alpha statistics and target other-background statistics. Before/after
covariance includes all shared inputs; statistical only, first-order propagation.
Three unit tests pass. Full archive/plot validation passes. Largest SF shift
in common chi2 bins is 0.365615 of the baseline prediction's statistical sigma.

Combined diagnostics sum separate-channel chi2s and yields without pooling alpha
or rebinning channels. Chi2 uses shared masks, target Data>=5 and DY Neff>=5 in
both regions/weight choices. No fitted parameters. Sparse bins remain plotted.

| Category | SR chi2/ndof noSF -> SF | VR chi2/ndof noSF -> SF |
| --- | ---: | ---: |
| Resolved untag | 130.74/62 -> same | 191.07/62 -> same |
| Resolved btag | 18.00/28 -> same | 34.31/33 -> same |
| Resolved VBF | 55.26/52 -> same | 54.27/55 -> same |
| Merged untag | 34.77/37 ->35.26/37 | 13.24/14 ->13.32/14 |
| Merged btag | 32.43/27 ->32.49/27 | N/A, sparse |
| Merged VBF | 16.77/18 ->16.93/18 | N/A, sparse |

No meaningful improvement in alpha closure from PNet alone; merged changes are
small compared with current statistics. Resolved untag nonclosure remains,
particularly VR: Data/Bkg=0.92184, and statistical yield pull=-12.09 in covered,
defined bins. SR resolved untag integrated ratio=0.99751 despite shape chi2.

Merged SR Data/Bkg noSF->SF: untag0.98233->0.98259, btag0.99157->0.99074,
VBF0.89468->0.89613. VR: untag1.15388->1.16024, btag1.21697->1.19879,
VBF1.40310->1.40992. These yields sum only covered bins with defined alpha;
undefined-alpha bins exclude 4, 8 and 1 observed merged VR events respectively.
They are not a complete full-search normalization test. Binwise Gaussian chi2
for merged VR btag/VBF is unavailable at this original grid: 0/0 in CSV is an
empty sum, not perfect closure. Human summary explicitly labels N/A.

Saved 24 individual PNG/PDF plots, four channel/region six-category overview
pairs, two twelve-page bundles, 72 compatibility rows, primitive/prediction
archives, metadata with input archive SHA256, documentation, tests and validator.
Ratio panels correctly label Data/Bkg rather than Data/MC for this data-derived
estimate. Cached replotting changes display only. Original analyses preserved.
