"""Shared helpers for the AlphaRatio examination scripts.

Extracted from FitTemplate/AlphaRatio.ipynb ("2D distributions" and
"Study different CRs" sections). Behaviour mirrors the notebook; the only
change is that the input/output directories are anchored to this file's
location so the scripts can be run from any working directory.
"""

import array
import os
from pathlib import Path

# lxlogin exposes many CPUs but enforces a low per-user process/thread limit.
for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np
import ROOT as rt

SCRIPT_DIR = Path(__file__).resolve().parent  # NewScheme/studies/alpha_ratio
FIT_DIR = Path(os.environ.get('HZZ_ANALYSIS_ROOT',SCRIPT_DIR.parents[2])) / "FitTemplate"
RED_TREES = FIT_DIR / "RedTrees"
RED_TREES_OLD = FIT_DIR / "RedTrees_old"
OUT_DIR = FIT_DIR / "AlphaRatio_examination" / "plots"

YEAR = 2018


def safeDivide(hnum, hden, name):
    hratio = hnum.Clone(name)
    hratio.Reset("ICES")
    hratio.SetDirectory(0)

    nbins = hnum.GetNbinsX()
    for i in range(1, nbins + 1):
        num = hnum.GetBinContent(i)
        den = hden.GetBinContent(i)
        num_err = hnum.GetBinError(i)
        den_err = hden.GetBinError(i)

        if num == 0 or den == 0:
            hratio.SetBinContent(i, 0.0)
            hratio.SetBinError(i, 0.0)
        else:
            val = num / den
            err = ((num_err / den) ** 2 + (num * den_err / (den * den)) ** 2) ** 0.5
            hratio.SetBinContent(i, val)
            hratio.SetBinError(i, err)

    return hratio


SIGMA_DYNAMIC_BINNING = {
    "2016": {
        "2e": {
            "resolved": [0, 250, 260, 270, 280, 290, 300, 314.83968276, 330.528212014, 347.072083828, 364.475240167, 382.739270141, 401.863642894, 421.845963693, 442.682245095, 464.367185796, 486.894450649, 510.256946358, 534.447088401, 559.457055803, 585.279031336, 611.905425588, 639.329084156, 667.54347781, 696.542876047, 726.322504838, 756.878689707, 788.208985519, 820.312294469, 853.188973913, 886.840935699, 921.271738672, 956.486676032, 992.492859172, 1184.77898268, 1399.10585347, 1588.71301428, 2028.5089526, 3500],
            "merged": [0, 300, 323.758415279, 348.733471772, 361.694744885, 374.981433873, 388.600886316, 402.56057734, 416.86811042, 431.531218128, 446.5577628, 461.955737138, 477.733264735, 493.898600519, 510.460131118, 527.426375134, 544.805983336, 562.607738754, 580.84055668, 599.513484565, 618.63570182, 638.216519502, 658.265379892, 678.79185596, 699.805650706, 721.316596382, 743.334653586, 765.869910225, 788.932580343, 812.533002813, 836.681639879, 861.389075562, 886.666013905, 912.523277071, 938.97180328, 966.022644582, 993.686964467, 1205.47447583, 1414.41383441, 1609.8620194, 2015.7933192, 3500]
        },
        "2mu": {
            "resolved": [0, 250, 260, 270, 280, 290, 300, 315.766185777, 332.080421856, 348.96175908, 366.429910521, 384.5052745, 403.208958406, 422.562803343, 442.589409639, 463.312163233, 484.755262986, 506.943748939, 529.903531552, 553.661421961, 578.24516329, 603.683463042, 630.006026627, 657.243592048, 685.427965794, 714.592059985, 744.769930802, 775.996818259, 808.30918735, 841.744770635, 876.342612298, 912.143113745, 949.188080779, 987.520772418, 1200.12240523, 1398.38781908, 1625.69197308, 2030.56748665, 3500],
            "merged": [0, 300, 323.018552543, 347.455306709, 360.232608265, 373.397639662, 386.962166642, 400.938311981, 415.338566323, 430.175799342, 445.463271246, 461.214644626, 477.443996669, 494.165831735, 511.395094322, 529.147182423, 547.437961285, 566.283777591, 585.701474066, 605.708404531, 626.322449412, 647.562031721, 669.44613352, 691.994312885, 715.226721387, 739.1641221, 763.827908156, 789.240121869, 815.423474429, 842.401366203, 870.197907647, 898.837940854, 928.347061754, 958.751642991, 990.078857495, 1199.04624, 1404.22131742, 1591.93308247, 1979.60833501, 3500]
        }
    },
    "2017": {
        "2e": {
            "resolved": [0, 250, 260, 270, 280, 290, 300, 315.972662605, 332.764508945, 350.38460801, 368.840839679, 388.140031422, 408.288102422, 429.290212187, 451.150911078, 473.874290584, 497.464131593, 521.924049255, 547.257633419, 573.468583925, 600.560840346, 628.538705971, 657.406966085, 687.171000727, 717.836892283, 749.411528353, 781.902700453, 815.319199146, 849.670906278, 884.968885021, 921.225468438, 958.454347347, 996.670658228, 1203.20601429, 1388.38345515, 1593.44289801, 2008.76153096, 3500],
            "merged": [0, 300, 325.748044588, 353.091512848, 367.385162305, 382.106314331, 397.26481655, 412.870657487, 428.933964876, 445.465003776, 462.474174499, 479.972010343, 497.969175121, 516.476460479, 535.504783006, 555.065181118, 575.168811724, 595.826946659, 617.05096888, 638.852368429, 661.242738149, 684.233769153, 707.837246036, 732.065041845, 756.929112775, 782.441492612, 808.614286909, 835.459666898, 862.989863126, 891.217158829, 920.153883025, 949.812403348, 980.205118593, 1011.344451, 1214.55593456, 1406.70472747, 1621.20691966, 2014.49918213, 3500]
        },
        "2mu": {
            "resolved": [0, 250, 260, 270, 280, 290, 300, 316.26200357, 333.110628696, 350.56703666, 368.653152098, 387.391690538, 406.806186926, 426.921025189, 447.761468858, 469.353692801, 491.724816094, 514.902936082, 538.917163672, 563.797659891, 589.575673771, 616.283581593, 643.954927552, 672.62446589, 702.328204541, 733.103450358, 764.988855973, 798.024468335, 832.251779016, 867.713776319, 904.454999268, 942.52159355, 981.961369474, 1201.55560873, 1407.50635101, 1582.30257444, 1992.91510505, 3500],
            "merged": [0, 300, 325.236246118, 352.089519863, 366.154666371, 380.663433251, 395.629812468, 411.068237299, 426.993596253, 443.421247427, 460.367033321, 477.84729611, 495.878893412, 514.479214538, 533.666197266, 553.458345139, 573.874745309, 594.935086943, 616.659680216, 639.069475889, 662.186085522, 686.031802312, 710.629622591, 736.003268006, 762.177208392, 789.176685375, 817.027736711, 845.757221397, 875.392845574, 905.963189246, 937.497733841, 970.026890643, 1003.58203012, 1187.91647008, 1403.21333012, 1601.21040713, 2012.78323841, 3500]
        }
    },
    "2018": {
        "2e": {
            "resolved": [0, 250, 260, 270, 280, 290, 300, 315.14574398, 331.219361909, 348.230469516, 366.185628321, 385.088581765, 404.940538335, 425.740490065, 447.485555062, 470.171333521, 493.792267907, 518.341999368, 543.813714021, 570.200474229, 597.495531465, 625.692618651, 654.786221, 684.771825338, 715.646148692, 747.407347486, 780.055209217, 813.591328739, 848.01927158, 883.344726795, 919.575651994, 956.72241316, 994.797921919, 1199.75440892, 1382.84756921, 1585.73842068, 2000.92316448, 3500],
            "merged": [0, 300, 324.956643936, 351.214899613, 364.852237629, 378.839381716, 393.184755285, 407.89692156, 422.984583933, 438.456586161, 454.321912393, 470.58968703, 487.269174411, 504.369778313, 521.901041281, 539.872643757, 558.294403022, 577.176271936, 596.528337479, 616.360819075, 636.684066703, 657.508558793, 678.844899882, 700.703818052, 723.096162117, 746.032898576, 769.525108309, 793.583983025, 818.220821442, 843.447025209, 869.274094554, 895.713623657, 922.777295746, 950.476877906, 978.824215604, 1007.83122692, 1196.40141198, 1411.77831749, 1613.64121428, 1982.88456439, 3500]
        },
        "2mu": {
            "resolved": [0, 250, 260, 270, 280, 290, 300, 316.436622949, 333.452823644, 351.069038785, 369.306425703, 388.186887759, 407.733100662, 427.968539693, 448.917507907, 470.605165316, 493.057559107, 516.30165493, 540.365369278, 565.277603018, 591.068276101, 617.768363496, 645.409932394, 674.026180715, 703.651476986, 734.321401613, 766.072789619, 798.943774876, 832.97383591, 868.203843316, 904.676108837, 942.434436192, 981.524173672, 1198.64635003, 1401.56581176, 1573.32022664, 1975.34184205, 3500],
            "merged": [0, 300, 324.755566579, 350.966677813, 364.644613617, 378.718914869, 393.201067627, 408.102890798, 423.436545783, 439.214546399, 455.449769097, 472.155463467, 489.34526305, 507.033196469, 525.233698871, 543.961623714, 563.232254885, 583.061319176, 603.464999116, 624.459946177, 646.063294369, 668.292674216, 691.166227148, 714.702620305, 738.921061775, 763.841316262, 789.483721223, 815.869203465, 843.019296219, 870.956156715, 899.702584268, 929.28203888, 959.718660388, 991.037288165, 1198.95556659, 1401.62283074, 1585.94135907, 2024.21474827, 3500]
        }
    }
}

for _year in SIGMA_DYNAMIC_BINNING:
    for _lep in SIGMA_DYNAMIC_BINNING[_year]:
        for _cat in SIGMA_DYNAMIC_BINNING[_year][_lep]:
            SIGMA_DYNAMIC_BINNING[_year][_lep][_cat] = array.array(
                "d", SIGMA_DYNAMIC_BINNING[_year][_lep][_cat]
            )


def make_hist(name, year, lep, cat):
    bins = SIGMA_DYNAMIC_BINNING[str(year)][lep][cat]
    h = rt.TH1F(name, name, len(bins) - 1, bins)
    h.SetDirectory(0)
    h.Sumw2()
    return h


_TREE_ARRAY_CACHE = {}


def clear_tree_array_cache():
    """Release arrays after one category/channel/tag combination is complete."""
    _TREE_ARRAY_CACHE.clear()


def tree_arrays(tree):
    """Read the three analysis branches once per TChain.

    Iterating over a TTree from Python is the dominant cost of these scripts.
    RDataFrame performs the read in C++ and returns NumPy arrays; the small cache
    then lets all CR windows reuse the same data.
    """
    key = id(tree)
    if key not in _TREE_ARRAY_CACHE:
        arrays = rt.RDataFrame(tree).AsNumpy(["mass", "mass2q", "weight"])
        _TREE_ARRAY_CACHE[key] = tuple(
            np.asarray(arrays[name], dtype=np.float64)
            for name in ("mass", "mass2q", "weight")
        )
    return _TREE_ARRAY_CACHE[key]


def _fill_arrays(h, masses, weights, reset=True):
    """Vectorized TH1 fill, including under/overflow and Sumw2 errors."""
    if reset:
        h.Reset("ICES")
    axis = h.GetXaxis()
    edges = np.fromiter(
        (axis.GetBinLowEdge(i) for i in range(1, h.GetNbinsX() + 2)),
        dtype=np.float64,
        count=h.GetNbinsX() + 1,
    )
    root_bins = np.searchsorted(edges, masses, side="right")
    root_bins[masses < edges[0]] = 0
    sums = np.bincount(root_bins, weights=weights, minlength=h.GetNbinsX() + 2)
    sumw2 = np.bincount(
        root_bins, weights=weights * weights, minlength=h.GetNbinsX() + 2
    )
    for bin_index in range(h.GetNbinsX() + 2):
        if not reset:
            sums[bin_index] += h.GetBinContent(bin_index)
            sumw2[bin_index] += h.GetBinError(bin_index) ** 2
        h.SetBinContent(bin_index, sums[bin_index])
        h.SetBinError(bin_index, np.sqrt(sumw2[bin_index]))


def outside_window(mass2q, cut_low, cut_high):
    """Match the notebook selection: reject only cut_low < mass2q < cut_high."""
    return (mass2q <= cut_low) | (mass2q >= cut_high)


def fill_hist(tree, h, cut_low, cut_high):
    masses, mass2q, weights = tree_arrays(tree)
    selected = outside_window(mass2q, cut_low, cut_high)
    _fill_arrays(h, masses[selected], weights[selected])


def fill_hist_alpharatio(data, tt, vz, h, alpha, cut_low=0, cut_high=1000):
    h.Reset("ICES")
    axis = alpha.GetXaxis()
    edges = np.fromiter(
        (axis.GetBinLowEdge(i) for i in range(1, alpha.GetNbinsX() + 2)),
        dtype=np.float64,
        count=alpha.GetNbinsX() + 1,
    )
    alpha_values = np.fromiter(
        (alpha.GetBinContent(i) for i in range(alpha.GetNbinsX() + 2)),
        dtype=np.float64,
        count=alpha.GetNbinsX() + 2,
    )

    first = True
    for tree, sign in ((data, 1.0), (tt, -1.0), (vz, -1.0)):
        masses, mass2q, weights = tree_arrays(tree)
        selected = outside_window(mass2q, cut_low, cut_high)
        selected_masses = masses[selected]
        alpha_bins = np.searchsorted(edges, selected_masses, side="right")
        alpha_bins[selected_masses < edges[0]] = 0
        _fill_arrays(
            h,
            selected_masses,
            sign * weights[selected] * alpha_values[alpha_bins],
            reset=first,
        )
        first = False
