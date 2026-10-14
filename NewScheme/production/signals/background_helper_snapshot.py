"""Produce grouped reduced background and data trees for the HZZ analysis."""

import argparse
import logging
import os
import re
from collections import defaultdict
from pathlib import Path

import awkward as ak
import numpy as np
import uproot
import yaml


ERAS = ["2016preAPV", "2016postAPV", "2017", "2018"]
REGIONS = ["CR", "SR", "VR"]
CHANNELS = {"ee": "isEE", "mumu": "isMuMu"}
CATEGORIES = ["resolved", "merged"]
TAGS = ["btag", "untag", "vbftag"]
PARTICLENET_SAMPLES = {
    "TTTo2L2Nu",
    "WWTo2L2Nu",
    "WZTo2Q2L",
    "ZZTo2Q2L",
}

BRANCHES = {
    "resolved": {
        "mass": "mass2l2jet",
        "m2l2q": "mass2l2jet",
        "pt2l": "pt2l",
        "mass2l": "mass2l",
        "pt2q": "pt2jet",
        "mass2q": "mass2jet",
        "kd": "KD_Zjj",
        "kd_jjvbf": "KD_jjVBF",
        "kd_jvbf": "KD_jVBF",
        "particleNetZvsQCD": "particleNetZvsQCD",
        "particleNetZbbvslight": "particleNetZbbvslight",
    },
    "merged": {
        "mass": "mass2lj",
        "m2l2q": "mass2lj",
        "pt2l": "pt2l",
        "mass2l": "mass2l",
        "pt2q": "ptmerged",
        "mass2q": "massmerged",
        "kd": "KD_ZJ",
        "kd_jjvbf": "KD_JJVBF",
        "kd_jvbf": "KD_JVBF",
        "particleNetZvsQCD": "particleNetZvsQCD",
        "particleNetZbbvslight": "particleNetZbbvslight",
    },
}

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG_DIR = SCRIPT_DIR / "cards"
DEFAULT_SF_DIR = SCRIPT_DIR / "sf"
LOGGER = logging.getLogger(__name__)


def needs_particle_net(sample):
    return sample.startswith("DY") or sample in PARTICLENET_SAMPLES


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--era", choices=ERAS + ["all"], default="all")
    parser.add_argument("--region", choices=REGIONS + ["all"], default="all")
    parser.add_argument(
        "--channel", choices=list(CHANNELS) + ["all"], default="all"
    )
    parser.add_argument("--cat", choices=CATEGORIES + ["all"], default="all")
    parser.add_argument("--tag", choices=TAGS + ["all"], default="all")
    parser.add_argument("--output-dir", default="reduced_trees")
    parser.add_argument(
        "--alpha-sr",
        default=None,
        help=(
            "SR alpha ROOT file or pattern. Patterns may contain {era} and "
            "{year}; {year} is 2016 for both APV eras."
        ),
    )
    parser.add_argument(
        "--alpha-vr",
        default=None,
        help="VR alpha ROOT file or pattern; supports {era} and {year}",
    )
    parser.add_argument("--config-dir", type=Path, default=DEFAULT_CONFIG_DIR)
    parser.add_argument("--sf-dir", type=Path, default=DEFAULT_SF_DIR)
    parser.add_argument(
        "--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO"
    )
    return parser.parse_args()


def selected(value, allowed):
    return allowed if value == "all" else [value]


def sample_group(sample):
    if sample == "Data":
        return "Data"
    if sample.startswith("DY"):
        return "DY"
    if sample in {"TTTo2L2Nu", "WWTo2L2Nu", "tZq"}:
        return "TTplusWW"
    if sample in {"WZTo2Q2L", "ZZTo2Q2L"}:
        return "VZ"
    raise ValueError(
        f"Sample {sample!r} in Samples_lists has no requested output group"
    )


def input_path(config, era, sample):
    if sample == "Data":
        return Path(config["ori_dir"]) / era / "Data" / "skimed" / (
            f"Data{era}UL_noDuplicates.root"
        )
    return (
        Path(config["ori_dir"])
        / era
        / config["samples_inf"][sample][0]
        / "skimed"
        / f"{sample}.root"
    )


def selection_branches(config, regions, channels, category, tags):
    names = set()
    for region in regions:
        for channel in channels:
            for tag in tags:
                expression = config["cut"][region][CHANNELS[channel]][category][tag]
                names.update(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expression))
    names.difference_update({"True", "False"})
    return names


def read_sample(config, era, sample, category, regions, channels, tags):
    path = input_path(config, era, sample)
    if not path.is_file():
        raise FileNotFoundError(f"Missing input file: {path}")

    group = sample_group(sample)
    needed = selection_branches(config, regions, channels, category, tags)
    needed.update(BRANCHES[category].values())
    if sample != "Data":
        needed.add("EventWeight")
        if category == "resolved" and group != "DY":
            needed.add("w_deepjet")
        if category == "merged" and needs_particle_net(sample):
            needed.update({"isbjet", "iscjet", "islightjet", "ptmerged"})

    with uproot.open(path) as root_file:
        if "passedEvents" not in root_file:
            raise KeyError(f"Missing passedEvents in {path}")
        tree = root_file["passedEvents"]
        available = set(tree.keys())
        missing = needed - available
        if missing:
            raise KeyError(f"Missing branches {sorted(missing)} in {path}")
        events = tree.arrays(sorted(needed), library="ak")

        sum_weight = None
        if sample != "Data":
            if "sumWeights" not in root_file:
                raise KeyError(f"Missing sumWeights in {path}")
            sum_weight = float(root_file["sumWeights"].to_boost().sum())
            if not np.isfinite(sum_weight) or sum_weight == 0.0:
                raise ValueError(f"Invalid sumWeights={sum_weight} in {path}")

    LOGGER.info("Read %s (%s): %d events", sample, group, len(events))
    return events, sum_weight


def particle_net_sf(era, events, tagger, scale_factors):
    if tagger == "ZvsQCD":
        field = "particleNetZvsQCD"
        b_wps = (
            {"LP": (0.90, 0.94), "MP": (0.94, 0.98), "HP": (0.98, 1.0)}
            if "2016" in era
            else {"LP": (0.90, 0.94), "MP": (0.94, 0.96), "HP": (0.96, 1.0)}
        )
        c_wps = {
            "LP": (0.90, 0.94),
            "MP": (0.94, 0.98),
            "HP": (0.98, 1.0),
        }
        light_wps = (
            {
                "LP": (0.90, 0.94),
                "MP": (0.94, 0.98),
                "HP": (0.98, 1.0),
            }
            if era == "2016postAPV"
            else {"LP": (0.90, 0.94), "HP": (0.94, 1.0)}
        )
    else:
        field = "particleNetZbbvslight"
        b_wps = c_wps = light_wps = {"LP": (0.0, 0.9), "HP": (0.9, 1.0)}

    result = np.zeros(len(events), dtype=np.float64)
    assigned = np.zeros(len(events), dtype=np.int32)
    flavor_wps = {
        "bhadron": ("isbjet", b_wps),
        "chadron": ("iscjet", c_wps),
        "lighthadron": ("islightjet", light_wps),
    }
    pt_ranges = {
        "200To400": (200.0, 400.0),
        "400To600": (400.0, 600.0),
        "600ToInf": (600.0, 9999999.0),
    }

    for flavor, (flavor_field, working_points) in flavor_wps.items():
        for purity, (low, high) in working_points.items():
            for pt_name, (pt_low, pt_high) in pt_ranges.items():
                mask = (
                    (events[field] > low)
                    & (events[field] <= high)
                    & (events["ptmerged"] >= pt_low)
                    & (events["ptmerged"] < pt_high)
                    & (events[flavor_field] == True)
                )
                selected_mask = ak.to_numpy(mask)
                assigned[selected_mask] += 1
                result[selected_mask] = float(
                    scale_factors[tagger][flavor][purity][pt_name]
                )
    if np.any(assigned != 1):
        raise ValueError(f"{era}/{tagger}: {np.count_nonzero(assigned != 1)} selected jets outside or ambiguously inside ParticleNet calibration coverage")
    if np.any(~np.isfinite(result)) or np.any(result <= 0):
        raise ValueError("Invalid ParticleNet scale factor")
    return result


def mc_weights(events, config, era, sample, group, category, tag, sum_weight, sfs,
               return_components=False):
    weights = ak.to_numpy(
        events["EventWeight"]
        * config["lumi"][era]
        * 1000.0
        * config["samples_inf"][sample][1]
        / sum_weight
    ).astype(np.float64, copy=True)

    qcd = np.ones(len(events), dtype=np.float64)
    light = qcd.copy()
    if category == "merged" and needs_particle_net(sample):
        qcd = particle_net_sf(era, events, "ZvsQCD", sfs)
        if tag in {"btag", "untag"}:
            light = particle_net_sf(era, events, "ZvsLight", sfs)
    elif category == "resolved" and group != "DY" and tag in {"btag", "untag"}:
        deepjet = ak.to_numpy(events["w_deepjet"])
        if np.any(deepjet == -999) or np.any(~np.isfinite(deepjet)):
            raise ValueError("Invalid DeepJet weight in selected resolved non-VBF events")
        weights *= deepjet
    pnet = qcd * light
    full = weights * pnet
    if np.any(~np.isfinite(full)):
        raise ValueError("Nonfinite event weight")
    if return_components:
        return dict(weight=full, weight_no_pnet=weights,
                    pnet_sf=pnet, pnet_sf_zvsqcd=qcd, pnet_sf_zvslight=light)
    return full


def reduced_record(events, category, weights):
    record = {
        output: ak.to_numpy(events[source])
        for output, source in BRANCHES[category].items()
    }
    if isinstance(weights, dict):
        record.update({key: np.asarray(value, dtype=np.float64)
                       for key, value in weights.items()})
    else:
        record["weight"] = np.asarray(weights, dtype=np.float64)
    return record


def concatenate(records):
    return {
        branch: np.concatenate([record[branch] for record in records])
        for branch in records[0]
    }


def output_name(process, channel, category, tag, era, region):
    return f"{process}_{channel}_{category}_{tag}_{era}_{region}.root"


def write_tree(output_dir, process, channel, category, tag, era, region, record):
    path = output_dir / output_name(process, channel, category, tag, era, region)
    with uproot.recreate(path) as output_file:
        output_file["tree"] = record
    LOGGER.info(
        "Wrote %s: entries=%d sum(weight)=%.8g",
        path,
        len(record["weight"]),
        np.sum(record["weight"], dtype=np.float64),
    )


def alpha_path(pattern, era):
    year = "2016" if era.startswith("2016") else era
    return Path(pattern.format(era=era, year=year))


def apply_alpha(record, histogram):
    edges = np.asarray(histogram.axes[0].edges, dtype=np.float64)
    values = np.asarray(histogram.values(), dtype=np.float64)
    indices = np.searchsorted(edges, record["mass"], side="right") - 1
    valid = (indices >= 0) & (indices < len(values))
    result = {branch: values_.copy() for branch, values_ in record.items()}
    result["weight"] = np.zeros(len(record["weight"]), dtype=np.float64)
    result["weight"][valid] = (
        record["weight"][valid] * values[indices[valid]]
    )
    return result, np.count_nonzero(~valid)


def make_zjet(store, alpha_pattern, era, category, region, channels, tags):
    path = alpha_path(alpha_pattern, era)
    if not path.is_file():
        raise FileNotFoundError(f"Missing {region} alpha file: {path}")

    with uproot.open(path) as alpha_file:
        for channel in channels:
            for tag in tags:
                parts = []
                for process, sign in (("Data", 1.0), ("TTplusWW", -1.0), ("VZ", -1.0)):
                    key = (process, channel, tag, "CR")
                    if key not in store:
                        raise RuntimeError(
                            f"Cannot build Zjet: {process} is absent for "
                            f"{era}/{category}/{channel}/{tag}/CR. Check Samples_lists."
                        )
                    part = concatenate(store[key])
                    part["weight"] = sign * part["weight"]
                    parts.append(part)

                signed_cr = concatenate(parts)
                prefix = "resolved" if category == "resolved" else "merged_tag"
                histogram_name = f"{prefix}_{CHANNELS[channel]}_{tag}"
                if histogram_name not in alpha_file:
                    raise KeyError(f"Missing {histogram_name} in {path}")
                histogram = alpha_file[histogram_name].to_hist()
                zjet, outside = apply_alpha(signed_cr, histogram)
                LOGGER.info(
                    "Applied nominal %s alpha from %s (%d events outside range)",
                    region,
                    path,
                    outside,
                )
                yield channel, tag, zjet


def process_category(config, era, category, regions, channels, tags, args):
    alpha_patterns = {"SR": args.alpha_sr, "VR": args.alpha_vr}
    zjet_regions = [region for region in regions if alpha_patterns.get(region)]
    read_regions = list(regions)
    if zjet_regions and "CR" not in read_regions:
        read_regions.append("CR")

    enabled = config.get("Samples_lists", [])
    if not enabled:
        raise RuntimeError("Samples_lists is empty")
    groups = {sample: sample_group(sample) for sample in enabled}
    LOGGER.info("Samples_lists groups for %s: %s", era, groups)
    if zjet_regions:
        missing = {"Data", "TTplusWW", "VZ"} - set(groups.values())
        if missing:
            raise RuntimeError(
                "Zjet was requested but Samples_lists is missing groups: "
                + ", ".join(sorted(missing))
            )

    need_sfs = category == "merged" and any(
        needs_particle_net(sample) for sample in enabled
    )
    sfs = None
    if need_sfs:
        sf_path = args.sf_dir / f"NetSF_signal_{era}.yml"
        with sf_path.open(encoding="utf-8") as stream:
            sfs = yaml.safe_load(stream)

    store = defaultdict(list)
    for sample in enabled:
        group = groups[sample]
        events, sum_weight = read_sample(
            config, era, sample, category, read_regions, channels, tags
        )
        for region in read_regions:
            for channel in channels:
                for tag in tags:
                    expression = config["cut"][region][CHANNELS[channel]][category][tag]
                    chosen = events[ak.numexpr.evaluate(expression, events)]
                    weights = (
                        {key: np.ones(len(chosen), dtype=np.float64)
                         for key in ("weight", "weight_no_pnet", "pnet_sf",
                                     "pnet_sf_zvsqcd", "pnet_sf_zvslight")}
                        if sample == "Data"
                        else mc_weights(
                            chosen,
                            config,
                            era,
                            sample,
                            group,
                            category,
                            tag,
                            sum_weight,
                            sfs,
                            return_components=True,
                        )
                    )
                    store[(group, channel, tag, region)].append(
                        reduced_record(chosen, category, weights)
                    )

    output_dir = Path(args.output_dir)
    for (process, channel, tag, region), records in store.items():
        if region in regions:
            write_tree(
                output_dir,
                process,
                channel,
                category,
                tag,
                era,
                region,
                concatenate(records),
            )

    for region in zjet_regions:
        for channel, tag, zjet in make_zjet(
            store,
            alpha_patterns[region],
            era,
            category,
            region,
            channels,
            tags,
        ):
            write_tree(
                output_dir, "Zjet", channel, category, tag, era, region, zjet
            )


def main():
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    eras = selected(args.era, ERAS)
    regions = selected(args.region, REGIONS)
    channels = selected(args.channel, list(CHANNELS))
    categories = selected(args.cat, CATEGORIES)
    tags = selected(args.tag, TAGS)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.alpha_sr is None and args.alpha_vr is None:
        LOGGER.info("No alpha files were supplied; Zjet outputs are omitted")

    for era in eras:
        config_path = args.config_dir / f"config_UL{era}.yml"
        LOGGER.info("Reading %s", config_path)
        with config_path.open(encoding="utf-8") as stream:
            config = yaml.safe_load(stream)
        for category in categories:
            process_category(
                config, era, category, regions, channels, tags, args
            )


if __name__ == "__main__":
    main()
