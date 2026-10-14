"""Mass spectra of DY vs (Data - TT - VZ) in bins of mass2q.

Replaces the former 2D colz diagnostic. For each category / channel / tag
(2018, CR+SR stacked) the mass distribution is shown in five mass2q slices
([0,70,105,135,165,inf]). Each slice gets a panel overlaying the simulated DY
(Z+jets) with the data-driven Z+jets estimate Data - TT/WW - VZ, plus their
ratio. Outputs go to plots/mass2q_slices/<cat>/.
"""

import argparse
import os

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_variable, "1")
import numpy as np
import ROOT as rt

from common import (
    OUT_DIR,
    RED_TREES_OLD,
    YEAR,
    _fill_arrays,
    clear_tree_array_cache,
    make_hist,
    safeDivide,
    tree_arrays,
)

CATS = ["resolved", "merged"]
LEPS = ["2e", "2mu"]
TAGKEYS = ["", "btag", "vbf"]
TAGS = {"": "untag", "btag": "btag", "vbf": "vbftag"}
CHANNELS = {"2e": "ee", "2mu": "mumu"}
MASS2Q_EDGES = [0.0, 70.0, 105.0, 135.0, 165.0, float("inf")]
MASS_RANGE = (300.0, 2000.0)
SLICE_LABELS = [
    "0 < M_{2q} < 70",
    "70 < M_{2q} < 105",
    "105 < M_{2q} < 135",
    "135 < M_{2q} < 165",
    "M_{2q} > 165",
]


def _build_chain(process, channel, cat, tag, year):
    chain = rt.TChain("tree", "tree")
    for region in ("CR", "SR"):
        path = RED_TREES_OLD / f"{process}_{channel}_{cat}_{tag}_{year}_{region}.root"
        if not path.is_file():
            raise FileNotFoundError(f"Missing input tree: {path}")
        if chain.Add(str(path)) != 1:
            raise OSError(f"ROOT could not add input tree: {path}")
    return chain


def _make_slice_hists(prefix, cat, lep, tag, year):
    hists = [
        make_hist(f"{prefix}_{year}_{cat}_{lep}_{tag}_{i}", year, lep, cat)
        for i in range(len(SLICE_LABELS))
    ]
    if cat == "merged":
        hists = [h.Rebin(2) for h in hists]
    return hists


def _fill_slices(tree, hists):
    masses, mass2q, weights = tree_arrays(tree)
    slice_indices = np.searchsorted(MASS2Q_EDGES, mass2q, side="right") - 1
    for idx, hist in enumerate(hists):
        selected = slice_indices == idx
        _fill_arrays(hist, masses[selected], weights[selected])


def _panel_max(xmin, xmax, *hists):
    ymax = 0.0
    for h in hists:
        axis = h.GetXaxis()
        for b in range(1, h.GetNbinsX() + 1):
            if xmin <= axis.GetBinCenter(b) <= xmax:
                ymax = max(ymax, h.GetBinContent(b) + h.GetBinError(b))
    return ymax if ymax > 0.0 else 1.0


def _plot_slices(dy, zjet, cat, lep, tag, year):
    n = len(SLICE_LABELS)
    rowh = 1.0 / n
    ratio_frac = 0.32

    canvas = rt.TCanvas(
        f"c_{year}_{cat}_{lep}_{tag}", f"c_{year}_{cat}_{lep}_{tag}", 1000, 1800
    )

    keep_alive = []
    for i in range(n):
        y_top = (n - i) / n
        y_bot = (n - 1 - i) / n
        split = y_bot + rowh * ratio_frac

        canvas.cd()
        main = rt.TPad(f"main_{i}", "", 0.0, split, 1.0, y_top)
        ratio = rt.TPad(f"ratio_{i}", "", 0.0, y_bot, 1.0, split)
        main.SetBottomMargin(0.02)
        main.SetTopMargin(0.03)
        main.SetLeftMargin(0.13)
        ratio.SetTopMargin(0.02)
        ratio.SetBottomMargin(0.55)
        ratio.SetLeftMargin(0.13)
        canvas.cd()
        main.Draw()
        ratio.Draw()
        keep_alive.extend([main, ratio])

        main.cd()
        zjet[i].SetLineColor(rt.kBlack)
        zjet[i].SetLineWidth(2)
        zjet[i].SetMarkerColor(rt.kBlack)
        zjet[i].SetMarkerStyle(20)
        zjet[i].SetMarkerSize(0.7)
        dy[i].SetLineColor(rt.kBlue)
        dy[i].SetLineWidth(2)

        frame = rt.TH1F(f"frame_{i}", "", 100, MASS_RANGE[0], MASS_RANGE[1])
        frame.SetStats(0)
        frame.SetTitle("")
        frame.GetXaxis().SetLabelSize(0)
        frame.GetXaxis().SetTitleSize(0)
        frame.GetYaxis().SetTitle("Events")
        frame.GetYaxis().SetTitleSize(0.09)
        frame.GetYaxis().SetLabelSize(0.07)
        frame.GetYaxis().SetTitleOffset(0.9)
        frame.SetMinimum(0.0)
        frame.SetMaximum(_panel_max(MASS_RANGE[0], MASS_RANGE[1], dy[i], zjet[i]) * 1.3)
        frame.Draw("AXIS")
        zjet[i].Draw("E SAME")
        dy[i].Draw("HIST SAME")

        label = rt.TLatex()
        label.SetNDC()
        label.SetTextSize(0.12)
        label.DrawLatex(0.17, 0.80, SLICE_LABELS[i])
        keep_alive.extend([frame, label])

        if i == 0:
            leg = rt.TLegend(0.52, 0.60, 0.92, 0.88)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.AddEntry(dy[i], "DY (Z+jets MC)", "l")
            leg.AddEntry(zjet[i], "Data - TT - VZ", "ep")
            leg.Draw()
            keep_alive.append(leg)

        ratio.cd()
        ratio_frame = rt.TH1F(f"rframe_{i}", "", 100, MASS_RANGE[0], MASS_RANGE[1])
        ratio_frame.SetStats(0)
        ratio_frame.SetTitle("")
        ratio_frame.GetYaxis().SetTitle("Ratio")
        ratio_frame.GetYaxis().SetTitleSize(0.14)
        ratio_frame.GetYaxis().SetTitleOffset(0.35)
        ratio_frame.GetYaxis().SetLabelSize(0.12)
        ratio_frame.GetYaxis().SetNdivisions(505)
        ratio_frame.GetXaxis().SetTitle("M_{2l2q} [GeV]")
        ratio_frame.GetXaxis().SetTitleSize(0.14)
        ratio_frame.GetXaxis().SetTitleOffset(1.1)
        ratio_frame.GetXaxis().SetLabelSize(0.12)
        ratio_frame.SetMinimum(0.5)
        ratio_frame.SetMaximum(1.5)
        ratio_frame.Draw("AXIS")

        h_ratio = safeDivide(zjet[i], dy[i], f"ratio_{year}_{cat}_{lep}_{tag}_{i}")
        h_ratio.SetLineColor(rt.kBlack)
        h_ratio.SetMarkerColor(rt.kBlack)
        h_ratio.SetMarkerStyle(20)
        h_ratio.SetMarkerSize(0.5)
        h_ratio.Draw("E SAME")
        keep_alive.extend([ratio_frame, h_ratio])

        line = rt.TLine(MASS_RANGE[0], 1.0, MASS_RANGE[1], 1.0)
        line.SetLineStyle(2)
        line.Draw("SAME")
        keep_alive.append(line)

    out_dir = OUT_DIR / "mass2q_slices" / cat
    out_dir.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(out_dir / f"Zjet_{year}_{cat}_{lep}_{TAGS[tag]}.png"))
    return canvas


def run(year=YEAR, cats=CATS, leps=LEPS, tagkeys=TAGKEYS):
    for cat in cats:
        for lep in leps:
            channel = CHANNELS[lep]
            for tag in tagkeys:
                tag_name = TAGS[tag]
                dy = _make_slice_hists("DY", cat, lep, tag_name, year)
                data = _make_slice_hists("Data", cat, lep, tag_name, year)
                tt = _make_slice_hists("TTplusWW", cat, lep, tag_name, year)
                vz = _make_slice_hists("VZ", cat, lep, tag_name, year)

                _fill_slices(_build_chain("DY", channel, cat, tag_name, year), dy)
                _fill_slices(_build_chain("Data", channel, cat, tag_name, year), data)
                _fill_slices(_build_chain("TTplusWW", channel, cat, tag_name, year), tt)
                _fill_slices(_build_chain("VZ", channel, cat, tag_name, year), vz)

                zjet = []
                for i in range(len(SLICE_LABELS)):
                    h = data[i].Clone(f"Zjet_{i}")
                    h.Add(tt[i], -1)
                    h.Add(vz[i], -1)
                    zjet.append(h)

                _plot_slices(dy, zjet, cat, lep, tag, year)
                print(f"{year} {cat} {lep} {tag_name}: done")
                clear_tree_array_cache()


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
    reverse_tags = {value: key for key, value in TAGS.items()}
    run(
        year=args.year,
        cats=args.category or CATS,
        leps=args.lepton or LEPS,
        tagkeys=[reverse_tags[tag] for tag in args.tag] if args.tag else TAGKEYS,
    )


if __name__ == "__main__":
    main()
