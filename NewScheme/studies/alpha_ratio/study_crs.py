"""Study alternative CR windows for the Z+jets alpha extrapolation.

Ports the "Study different CRs" section of FitTemplate/AlphaRatio.ipynb. For each
candidate excluded-mass2q window it builds

    alpha(mass) = DY(SR) / DY(CR outside window)

and then compares the resulting CR-based Z+jets estimate,
alpha * [Data(CR) - TT/WW(CR) - WZ/ZZ(CR)], with the SR data-driven Z+jets
estimate Data(SR) - TT/WW(SR) - WZ/ZZ(SR).

Note (preserved from the notebook): SR trees are read from RedTrees/ while CR
trees are read from RedTrees_old/.
"""

import argparse
import os

# Import NumPy before ROOT so PyROOT's import hook cannot intercept NumPy's
# extension-module initialization in the older unbin environment.
for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_variable, "1")
import numpy as np  # noqa: F401
import ROOT as rt

from common import (
    OUT_DIR,
    RED_TREES,
    RED_TREES_OLD,
    SIGMA_DYNAMIC_BINNING,
    YEAR,
    clear_tree_array_cache,
    fill_hist,
    fill_hist_alpharatio,
    make_hist,
    safeDivide,
)

CATS = ["resolved", "merged"]
LEPS = ["2e", "2mu"]
TAGKEYS = ["", "btag", "vbf"]
PROCESSES = ["DY", "TTplusWW", "VZ", "Data"]
CHANNELS = {"2e": "ee", "2mu": "mumu"}
TAGS = {"": "untag", "btag": "btag", "vbf": "vbftag"}
CUTS = [[135, 165], [105, 165]]


def build_trees(year, cats=CATS, leps=LEPS, tagkeys=TAGKEYS):
    """Chain SR trees from RedTrees/ and CR trees from RedTrees_old/."""
    trees = {}
    for reg in ["SR", "CR"]:
        base = RED_TREES if reg == "SR" else RED_TREES_OLD
        for cat in cats:
            for lep in leps:
                channel = CHANNELS[lep]
                for tag in tagkeys:
                    tag_name = TAGS[tag]
                    for process in PROCESSES:
                        key = (process, year, reg, cat, lep, tag)
                        path = base / f"{process}_{channel}_{cat}_{tag_name}_{year}_{reg}.root"
                        if not path.is_file():
                            raise FileNotFoundError(f"Missing input tree: {path}")
                        trees[key] = rt.TChain("tree", "tree")
                        if trees[key].Add(str(path)) != 1:
                            raise OSError(f"ROOT could not add input tree: {path}")
    return trees


def build_alpha_ratio(
    trees, year, cuts=CUTS, cats=CATS, leps=LEPS, tagkeys=TAGKEYS
):
    """alpha = DY(SR) / DY(CR with mass2q outside the window)."""
    alpha_ratio = {}
    for cat in cats:
        for lep in leps:
            lep_name = "isEE" if lep == "2e" else "isMuMu"
            bins = SIGMA_DYNAMIC_BINNING[str(year)][lep][cat]
            for tag in tagkeys:
                tag_name = TAGS[tag]
                for cut_low, cut_high in cuts:
                    name = f"{cat}_{lep_name}_{tag_name}_{cut_low}_{cut_high}"
                    hnum = make_hist(f"hnum_{name}", year, lep, cat)
                    hden = make_hist(f"hden_{name}", year, lep, cat)
                    if cat == "merged":
                        hnum = hnum.Rebin(2)
                        hden = hden.Rebin(2)

                    fill_hist(
                        trees["DY", year, "SR", cat, lep, tag], hnum, -1, -1
                    )
                    fill_hist(
                        trees["DY", year, "CR", cat, lep, tag],
                        hden,
                        cut_low,
                        cut_high,
                    )

                    hist_name = (
                        name
                        if cat == "resolved"
                        else f"{cat}_tag_{lep_name}_{tag_name}_{cut_low}_{cut_high}"
                    )
                    alpha_ratio[cat, lep, tag, cut_low, cut_high] = safeDivide(
                        hnum, hden, hist_name
                    )
                    print(
                        f"cat={cat}, lep={lep}, tag={tag}, cut_low={cut_low}, "
                        f"cut_high={cut_high}: num {hnum.Integral()}, den {hden.Integral()}"
                    )
    return alpha_ratio


def build_cr_histograms(
    trees, alpha_ratio, year, cuts=CUTS, cats=CATS, leps=LEPS, tagkeys=TAGKEYS
):
    """Fill CR Z+jets (alpha-weighted) and SR data-driven Z+jets histograms."""
    data_cr = {}
    tt_cr = {}
    vz_cr = {}
    zjet = {}
    data_sr = {}

    for cat in cats:
        for lep in leps:
            for tag in tagkeys:
                tag_name = TAGS[tag]
                for ic, (cut_low, cut_high) in enumerate(cuts):
                    key = (year, cat, lep, tag, ic)

                    data_cr[key] = make_hist(
                        f"Data_{year}_{cat}_{lep}_{tag_name}_{ic}", year, lep, cat
                    )
                    tt_cr[key] = make_hist(
                        f"TTplusWW_{year}_{cat}_{lep}_{tag_name}_{ic}", year, lep, cat
                    )
                    vz_cr[key] = make_hist(
                        f"VZ_{year}_{cat}_{lep}_{tag_name}_{ic}", year, lep, cat
                    )
                    zjet[key] = make_hist(
                        f"Zjet_{year}_{cat}_{lep}_{tag_name}_{ic}", year, lep, cat
                    )

                    if cat == "merged":
                        data_cr[key].Rebin(2)
                        tt_cr[key].Rebin(2)
                        vz_cr[key].Rebin(2)
                        zjet[key].Rebin(2)

                    fill_hist(
                        trees["Data", year, "CR", cat, lep, tag],
                        data_cr[key],
                        cut_low,
                        cut_high,
                    )
                    fill_hist(
                        trees["TTplusWW", year, "CR", cat, lep, tag],
                        tt_cr[key],
                        cut_low,
                        cut_high,
                    )
                    fill_hist(
                        trees["VZ", year, "CR", cat, lep, tag],
                        vz_cr[key],
                        cut_low,
                        cut_high,
                    )
                    fill_hist_alpharatio(
                        trees["Data", year, "CR", cat, lep, tag],
                        trees["TTplusWW", year, "CR", cat, lep, tag],
                        trees["VZ", year, "CR", cat, lep, tag],
                        zjet[key],
                        alpha_ratio[cat, lep, tag, cut_low, cut_high],
                        cut_low,
                        cut_high,
                    )

                    print(key, "Zjet old =", zjet[key].Integral())

    for cat in cats:
        for lep in leps:
            for tag in tagkeys:
                tag_name = TAGS[tag]
                key = (year, cat, lep, tag)

                data_sr[key] = make_hist(
                    f"Data_SR_{year}_{cat}_{lep}_{tag_name}", year, lep, cat
                )
                tt_sr = make_hist(
                    f"TTplusWW_SR_{year}_{cat}_{lep}_{tag_name}", year, lep, cat
                )
                vz_sr = make_hist(
                    f"VZ_SR_{year}_{cat}_{lep}_{tag_name}", year, lep, cat
                )

                if cat == "merged":
                    data_sr[key].Rebin(2)
                    tt_sr.Rebin(2)
                    vz_sr.Rebin(2)

                fill_hist(trees["Data", year, "SR", cat, lep, tag], data_sr[key], -1, -1)
                fill_hist(trees["TTplusWW", year, "SR", cat, lep, tag], tt_sr, -1, -1)
                fill_hist(trees["VZ", year, "SR", cat, lep, tag], vz_sr, -1, -1)
                data_sr[key].Add(tt_sr, -1)
                data_sr[key].Add(vz_sr, -1)

    return zjet, data_sr


def plot_overlay(
    zjet, data_sr, year, cuts=CUTS, cats=CATS, leps=LEPS, tagkeys=TAGKEYS
):
    for cat in cats:
        for lep in leps:
            for tag in tagkeys:
                canvas = rt.TCanvas(
                    f"c_{year}_{cat}_{lep}_{tag}",
                    f"c_{year}_{cat}_{lep}_{tag}",
                    800,
                    800,
                )
                colors = [rt.kRed, rt.kBlue]
                for ic in [1, 0]:
                    key = (year, cat, lep, tag, ic)
                    zjet[key].SetLineColor(colors[ic])
                    zjet[key].SetLineWidth(2)

                    if ic == 1:
                        zjet[key].SetTitle(";M_{2l2q} [GeV];Events")
                        zjet[key].Draw("E")
                    else:
                        zjet[key].Draw("E SAME")

                data = data_sr[year, cat, lep, tag]
                data.SetLineColor(rt.kBlack)
                data.SetMarkerStyle(20)
                data.SetMarkerSize(1)
                data.Draw("E SAME")

                leg = rt.TLegend(0.6, 0.7, 0.9, 0.9)
                for ic in range(2):
                    leg.AddEntry(
                        zjet[year, cat, lep, tag, ic],
                        "Z+jets, exclude [{0},{1}]".format(
                            cuts[ic][0], cuts[ic][1]
                        ),
                        "l",
                    )
                leg.AddEntry(data, "Data")
                leg.Draw()

                canvas.SaveAs(str(OUT_DIR / f"Zjet_{year}_{cat}_{lep}_{tag}.png"))


def plot_ratio_panel(
    zjet, data_sr, year, cuts=CUTS, cats=CATS, leps=LEPS, tagkeys=TAGKEYS
):
    for cat in cats:
        for lep in leps:
            for tag in tagkeys:
                canvas = rt.TCanvas(
                    f"c_{year}_{cat}_{lep}_{tag}",
                    f"c_{year}_{cat}_{lep}_{tag}",
                    800,
                    800,
                )

                suffix = f"{year}_{cat}_{lep}_{tag or 'untag'}"
                pad1 = rt.TPad(f"pad1_{suffix}", "pad1", 0.0, 0.30, 1.0, 1.0)
                pad2 = rt.TPad(f"pad2_{suffix}", "pad2", 0.0, 0.00, 1.0, 0.30)
                pad1.SetBottomMargin(0.02)
                pad2.SetTopMargin(0.03)
                pad2.SetBottomMargin(0.35)
                pad1.Draw()
                pad2.Draw()

                colors = [rt.kRed, rt.kBlue]
                data_key = (year, cat, lep, tag)

                pad1.cd()
                for ic in [1, 0]:
                    key = (year, cat, lep, tag, ic)
                    zjet[key].SetLineColor(colors[ic])
                    zjet[key].SetLineWidth(2)

                    if ic == 1:
                        zjet[key].SetTitle(";M_{2l2q} [GeV];Events")
                        zjet[key].GetXaxis().SetLabelSize(0)
                        zjet[key].GetXaxis().SetTitleSize(0)
                        zjet[key].Draw("E")
                    else:
                        zjet[key].Draw("E SAME")

                data_sr[data_key].SetLineColor(rt.kBlack)
                data_sr[data_key].SetMarkerStyle(20)
                data_sr[data_key].SetMarkerSize(1)
                data_sr[data_key].Draw("E SAME")

                leg = rt.TLegend(0.6, 0.7, 0.9, 0.9)
                for ic in range(2):
                    leg.AddEntry(
                        zjet[year, cat, lep, tag, ic],
                        "Z+jets, exclude [{0},{1}]".format(
                            cuts[ic][0], cuts[ic][1]
                        ),
                        "l",
                    )
                leg.AddEntry(data_sr[data_key], "Data", "lep")
                leg.Draw()

                pad2.cd()
                pad2.SetGridy()

                ratios = []
                for ic in [1, 0]:
                    key = (year, cat, lep, tag, ic)
                    ratio = safeDivide(
                        data_sr[data_key],
                        zjet[key],
                        f"ratio_{year}_{cat}_{lep}_{tag}_{ic}",
                    )
                    ratio.SetLineColor(colors[ic])
                    ratio.SetMarkerColor(colors[ic])
                    ratio.SetMarkerStyle(20)
                    ratio.SetMarkerSize(0.7)
                    ratio.SetTitle("")
                    ratio.GetYaxis().SetTitle("Data / Z+jets")
                    ratio.GetXaxis().SetTitle("M_{2l2q} [GeV]")
                    ratio.SetMinimum(0.5)
                    ratio.SetMaximum(1.5)
                    ratio.GetYaxis().SetNdivisions(505)
                    ratio.GetYaxis().SetTitleSize(0.10)
                    ratio.GetYaxis().SetTitleOffset(0.45)
                    ratio.GetYaxis().SetLabelSize(0.09)
                    ratio.GetXaxis().SetTitleSize(0.12)
                    ratio.GetXaxis().SetTitleOffset(1.05)
                    ratio.GetXaxis().SetLabelSize(0.10)

                    if ic == 1:
                        ratio.Draw("E")
                    else:
                        ratio.Draw("E SAME")
                    ratios.append(ratio)

                xmin = data_sr[data_key].GetXaxis().GetXmin()
                xmax = data_sr[data_key].GetXaxis().GetXmax()
                line = rt.TLine(xmin, 1.0, xmax, 1.0)
                line.SetLineStyle(2)
                line.Draw("SAME")

                canvas.SaveAs(
                    str(OUT_DIR / f"Zjet_{year}_{cat}_{lep}_{tag}_inSR.png")
                )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, default=YEAR)
    parser.add_argument("--category", choices=CATS, action="append")
    parser.add_argument("--lepton", choices=LEPS, action="append")
    parser.add_argument("--tag", choices=list(TAGS.values()), action="append")
    return parser.parse_args()


def main():
    args = parse_args()
    rt.gROOT.SetBatch(True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cats = args.category or CATS
    leps = args.lepton or LEPS
    reverse_tags = {value: key for key, value in TAGS.items()}
    tagkeys = [reverse_tags[tag] for tag in args.tag] if args.tag else TAGKEYS
    # Process one combination at a time so the array cache is fast but bounded.
    for cat in cats:
        for lep in leps:
            for tag in tagkeys:
                trees = build_trees(args.year, [cat], [lep], [tag])
                alpha_ratio = build_alpha_ratio(
                    trees,
                    args.year,
                    cats=[cat],
                    leps=[lep],
                    tagkeys=[tag],
                )
                zjet, data_sr = build_cr_histograms(
                    trees,
                    alpha_ratio,
                    args.year,
                    cats=[cat],
                    leps=[lep],
                    tagkeys=[tag],
                )
                plot_overlay(
                    zjet,
                    data_sr,
                    args.year,
                    cats=[cat],
                    leps=[lep],
                    tagkeys=[tag],
                )
                plot_ratio_panel(
                    zjet,
                    data_sr,
                    args.year,
                    cats=[cat],
                    leps=[lep],
                    tagkeys=[tag],
                )
                clear_tree_array_cache()


if __name__ == "__main__":
    main()
