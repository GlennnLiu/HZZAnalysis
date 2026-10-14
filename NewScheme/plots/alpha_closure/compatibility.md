# Alpha-method SR/VR compatibility, 2018

Combined tests sum separate ee/mumu chi2 values; original channel binning. Statistical only.

| Region | Category | chi2/ndof no PNet | with PNet | Data/Bkg no PNet | with PNet |
| --- | --- | ---: | ---: | ---: | ---: |
| SR | resolved untag | 130.74/62 | 130.74/62 | 0.9975 | 0.9975 |
| SR | resolved btag | 18.00/28 | 18.00/28 | 0.9811 | 0.9811 |
| SR | resolved vbftag | 55.26/52 | 55.26/52 | 1.0194 | 1.0194 |
| SR | merged untag | 34.77/37 | 35.26/37 | 0.9823 | 0.9826 |
| SR | merged btag | 32.43/27 | 32.49/27 | 0.9916 | 0.9907 |
| SR | merged vbftag | 16.77/18 | 16.93/18 | 0.8947 | 0.8961 |
| VR | resolved untag | 191.07/62 | 191.07/62 | 0.9218 | 0.9218 |
| VR | resolved btag | 34.31/33 | 34.31/33 | 0.8841 | 0.8841 |
| VR | resolved vbftag | 54.27/55 | 54.27/55 | 0.9107 | 0.9107 |
| VR | merged untag | 13.24/14 | 13.32/14 | 1.1539 | 1.1602 |
| VR | merged btag | N/A (sparse) | N/A (sparse) | 1.2170 | 1.1988 |
| VR | merged vbftag | N/A (sparse) | N/A (sparse) | 1.4031 | 1.4099 |

Yield ratios include only covered bins with defined alpha. Undefined-alpha Data yields are explicit in CSV;
in merged VR they exclude 4 untag, 8 btag and 1 VBF events. No Gaussian chi2 is available for
merged VR btag/VBF at original binning: CSV chi2=0, ndof=0 is an empty sum, not perfect closure.

Max SF shift / baseline prediction statistical sigma on chi2 bins: 0.3656.

No calibration/modeling systematics; alternative-method chi2 differences are not independent significance tests.
