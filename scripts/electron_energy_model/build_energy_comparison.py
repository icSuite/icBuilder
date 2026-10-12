"""Build a compact DMSP electron-energy model comparison dataset.

The detector crossing files contain one row per DMSP second and IMAGE frame.
This script reduces those rows to one physically meaningful DMSP energy per
IMAGE-frame/satellite window.  It then evaluates the grid-independent
Zhang--Paxton latitude collapse at the matched Kp and circular-mean MLT.

Every crossing is retained.  Manual annotation is stored as a status field so
the same output can support both an all-crossing analysis and an
accepted-crossing analysis without rebuilding the scientific inputs.
"""

#%% Imports

import argparse
from pathlib import Path

from icbuilder.kp import DEFAULT_KP_PATH, load_gfz_kp, match_gfz_kp
from icphysics import collapse_zhang_paxton, hardy_electron_precipitation
import numpy as np
import pandas as pd
from tqdm import tqdm
import xarray as xr
from zhangpaxton2008 import zhang_paxton


#%% Defaults and input contract

REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_MATCHES = Path(
    "/home/bing/Dropbox/work/temp_storage/icBuilder_pipeline_test/"
    "matched_dmsp_image_data"
)
DEFAULT_ANNOTATIONS = REPOSITORY / "data/dmsp_frame_annotations.csv"
DEFAULT_OUTPUT = REPOSITORY / "data/electron_energy_model/comparison_data.nc"

DEFAULT_SUPPORT_SECONDS = (66, 120)
MAX_ENERGY_FRACTIONAL_STD = 0.25
MAX_FLUX_FRACTIONAL_STD = 0.20
ZP_MLT_STEP_HOURS = 0.05

INPUT_COLUMNS = [
    "orbit", "image_frame", "image_time", "dmsp_sat", "dmsp_time",
    "time_offset_seconds", "dmsp_mlat", "dmsp_mlt",
    "detector_separation_deg", "dmsp_electron_mean_energy",
    "dmsp_electron_mean_energy_fractional_std",
    "dmsp_electron_total_energy_flux",
    "dmsp_electron_total_energy_flux_fractional_std",
    "img_ratio", "img_energy", "img_energy_std", "method_valid"
]

ZP_FIELDS = {
    "representative_energy": "zp_area_mean_energy",
    "area_weighted_median_energy": "zp_area_median_energy",
    "q_weighted_energy": "zp_flux_weighted_energy",
    "weighted_spread": "zp_latitude_spread",
    "selected_lower_mlat": "zp_selected_lower_mlat",
    "selected_upper_mlat": "zp_selected_upper_mlat",
    "selected_area_fraction": "zp_selected_area_fraction"
}


#%% Annotation keys

def load_annotation_status(filename):
    """Return manual status keyed by orbit, image time, and satellite."""

    annotations = pd.read_csv(filename)
    required = {"orbit", "img_time", "satellites", "accepted"}
    missing = required.difference(annotations.columns)
    if missing:
        raise ValueError(f"{filename} is missing {sorted(missing)}")

    annotations["image_time"] = pd.to_datetime(annotations["img_time"])
    annotations["dmsp_sat"] = annotations["satellites"].str.lower()
    keys = ["orbit", "image_time", "dmsp_sat"]
    if annotations.duplicated(keys).any():
        raise ValueError("annotation orbit/time/satellite keys must be unique")

    return {
        (int(row.orbit), np.datetime64(row.image_time, "ns"), row.dmsp_sat):
            int(row.accepted)
        for row in annotations.itertuples()
    }


#%% DMSP window reduction

def circular_mean_mlt(values):
    """Return the circular mean magnetic local time in hours."""

    angle = np.asarray(values, dtype=float) * 2 * np.pi / 24
    mean_angle = np.arctan2(np.nanmean(np.sin(angle)), np.nanmean(np.cos(angle)))
    return float(np.mod(mean_angle * 24 / (2 * np.pi), 24))


def reduce_crossing(group, support_seconds, annotation_status):
    """Reduce one IMAGE-frame/satellite crossing at one temporal support."""

    half_support = support_seconds / 2
    sample = group.loc[
        group["time_offset_seconds"].ge(-half_support)
        & group["time_offset_seconds"].lt(half_support)
    ]

    good = (
        np.isfinite(sample["dmsp_electron_mean_energy"])
        & sample["dmsp_electron_mean_energy"].gt(0)
        & np.isfinite(sample["dmsp_electron_total_energy_flux"])
        & sample["dmsp_electron_total_energy_flux"].gt(0)
        & np.isfinite(sample["dmsp_electron_mean_energy_fractional_std"])
        & sample["dmsp_electron_mean_energy_fractional_std"].le(
            MAX_ENERGY_FRACTIONAL_STD
        )
        & np.isfinite(
            sample["dmsp_electron_total_energy_flux_fractional_std"]
        )
        & sample["dmsp_electron_total_energy_flux_fractional_std"].le(
            MAX_FLUX_FRACTIONAL_STD
        )
    )
    dmsp = sample.loc[good]

    orbit = int(group["orbit"].iloc[0])
    image_frame = int(group["image_frame"].iloc[0])
    image_time = np.datetime64(group["image_time"].iloc[0], "ns")
    satellite = str(group["dmsp_sat"].iloc[0]).lower()
    key = (orbit, image_time, satellite)

    result = {
        "orbit": orbit,
        "image_frame": image_frame,
        "image_time": image_time,
        "year": int(str(image_time)[:4]),
        "satellite": satellite,
        "support_seconds": int(support_seconds),
        "annotation_status": int(annotation_status.get(key, -1)),
        "n_available": int(len(sample)),
        "n_dmsp": int(len(dmsp)),
        "dmsp_coverage": min(len(dmsp) / support_seconds, 1.0)
    }
    if len(dmsp) == 0:
        return result

    energy = dmsp["dmsp_electron_mean_energy"].to_numpy(dtype=float)
    flux = dmsp["dmsp_electron_total_energy_flux"].to_numpy(dtype=float)
    result.update({
        # Sum energy flux divided by the corresponding number-flux proxy.
        "dmsp_energy": float(np.sum(flux) / np.sum(flux / energy)),
        "dmsp_flux": float(np.mean(flux)),
        "dmsp_energy_fractional_std": float(np.median(
            dmsp["dmsp_electron_mean_energy_fractional_std"]
        )),
        "dmsp_flux_fractional_std": float(np.median(
            dmsp["dmsp_electron_total_energy_flux_fractional_std"]
        )),
        "mlt": circular_mean_mlt(dmsp["dmsp_mlt"]),
        "mlat": float(np.mean(dmsp["dmsp_mlat"])),
        "detector_separation_median_deg": float(np.median(
            dmsp["detector_separation_deg"]
        )),
        "detector_separation_maximum_deg": float(np.max(
            dmsp["detector_separation_deg"]
        ))
    })

    pointwise_good = (
        np.isfinite(dmsp["zp_pointwise_energy"])
        & dmsp["zp_pointwise_energy"].gt(0)
        & np.isfinite(dmsp["zp_pointwise_flux"])
        & dmsp["zp_pointwise_flux"].gt(0)
    )
    pointwise = dmsp.loc[pointwise_good]
    result["n_zp_pointwise"] = int(len(pointwise))
    if len(pointwise):
        pointwise_energy = pointwise["zp_pointwise_energy"].to_numpy(dtype=float)
        pointwise_flux = pointwise["zp_pointwise_flux"].to_numpy(dtype=float)
        result.update({
            # Apply the same energy-flux/number-flux reduction used for DMSP.
            "zp_pointwise_energy": float(
                np.sum(pointwise_flux)
                / np.sum(pointwise_flux / pointwise_energy)
            ),
            "zp_pointwise_flux": float(np.mean(pointwise_flux))
        })

    # Reduce independent Hardy moments over exactly the same DMSP window.
    if "hardy_pointwise_energy" in dmsp:
        hardy_good = np.isfinite(dmsp["hardy_pointwise_energy"]) & dmsp["hardy_pointwise_energy"].gt(0)
        hardy_good &= np.isfinite(dmsp["hardy_pointwise_flux"]) & dmsp["hardy_pointwise_flux"].gt(0)
        hardy = dmsp.loc[hardy_good]
        result["n_hardy_pointwise"] = len(hardy)
        result["hardy_pointwise_energy"] = np.nan
        result["hardy_pointwise_flux"] = np.nan
        if len(hardy):
            energy = hardy["hardy_pointwise_energy"].to_numpy(dtype=float)
            flux = hardy["hardy_pointwise_flux"].to_numpy(dtype=float)
            result["hardy_pointwise_energy"] = float(np.sum(flux) / np.sum(flux / energy))
            result["hardy_pointwise_flux"] = float(np.mean(flux))

    image_good = (
        dmsp["method_valid"].astype(bool)
        & np.isfinite(dmsp["img_energy"])
        & dmsp["img_energy"].gt(0)
    )
    image = dmsp.loc[image_good]
    result["n_image_energy"] = int(len(image))
    if len(image):
        # This is a descriptive along-track summary, not a fitted predictor.
        result.update({
            "image_ratio_energy_median": float(np.median(image["img_energy"])),
            "image_ratio_energy_std_median": float(np.median(
                image["img_energy_std"]
            )),
            "image_ratio_median": float(np.median(image["img_ratio"]))
        })

    return result


def reduce_orbit(filename, support_seconds, annotation_status, kp_series):
    """Read and reduce one detector-crossing orbit file."""

    with xr.open_dataset(filename) as data:
        if data.sizes.get("sample", 0) == 0:
            return []
        missing = set(INPUT_COLUMNS).difference(data.variables)
        if missing:
            raise ValueError(f"{filename} is missing {sorted(missing)}")

        available_half_support = float(
            data.attrs.get("frame_half_support_seconds", 0)
        )
        if max(support_seconds) / 2 > available_half_support:
            raise ValueError(
                f"{filename} contains only +/-{available_half_support:g} s, "
                f"not the requested {max(support_seconds):g}-s support"
            )
        matches = data[INPUT_COLUMNS].to_dataframe().reset_index(drop=True)

    # Evaluate the original model at every one-second DMSP position. This is
    # intentionally done before temporal reduction; evaluating once at the
    # mean latitude would not describe the same satellite window as DMSP.
    matched_kp = match_gfz_kp(matches["image_time"].to_numpy(), kp_series)
    pointwise_energy, pointwise_flux = zhang_paxton(
        matched_kp["kp"],
        matches["dmsp_mlt"].to_numpy(dtype=float),
        matches["dmsp_mlat"].to_numpy(dtype=float)
    )
    matches["zp_pointwise_energy"] = pointwise_energy
    matches["zp_pointwise_flux"] = pointwise_flux
    hardy = hardy_electron_precipitation(matched_kp["kp"], matches["dmsp_mlt"].to_numpy(dtype=float), matches["dmsp_mlat"].to_numpy(dtype=float))
    matches["hardy_pointwise_energy"] = hardy["mean_energy"]
    matches["hardy_pointwise_flux"] = hardy["energy_flux"]

    records = []
    keys = ["orbit", "image_frame", "image_time", "dmsp_sat"]
    groups = matches.groupby(keys, sort=False, observed=True)
    for _, group in groups:
        for seconds in support_seconds:
            records.append(reduce_crossing(group, seconds, annotation_status))
    return records


#%% Zhang--Paxton evaluation

def add_zhang_paxton(frame, kp_path, mlt_step=ZP_MLT_STEP_HOURS):
    """Attach Kp and periodic interpolations of all collapse energies."""

    kp_series = load_gfz_kp(kp_path)
    matched = match_gfz_kp(frame["image_time"].to_numpy(), kp_series)
    frame["kp"] = matched["kp"]
    frame["kp_interval_start"] = matched["interval_start"]
    for output_name in ZP_FIELDS.values():
        frame[output_name] = np.nan

    mlt_grid = np.arange(0, 24, mlt_step)
    extended_mlt = np.append(mlt_grid, 24)
    for kp in tqdm(np.sort(frame["kp"].unique()), desc="Collapse Zhang--Paxton"):
        collapsed = collapse_zhang_paxton(np.full(mlt_grid.shape, kp), mlt_grid)
        choose = frame["kp"].eq(kp) & np.isfinite(frame["mlt"])
        sample_mlt = np.mod(frame.loc[choose, "mlt"].to_numpy(), 24)

        for source_name, output_name in ZP_FIELDS.items():
            values = np.asarray(collapsed[source_name], dtype=float)
            periodic_values = np.append(values, values[0])
            frame.loc[choose, output_name] = np.interp(
                sample_mlt, extended_mlt, periodic_values
            )

    return frame


#%% Serialization and command line

def save_comparison(frame, filename, matches_path, annotations_path, kp_path,
                    support_seconds):
    """Write the compact comparison table with its scientific choices."""

    frame = frame.sort_values(
        ["image_time", "satellite", "support_seconds"]
    ).reset_index(drop=True)
    dataset = frame.to_xarray().rename({"index": "record"})
    dataset.attrs.update({
        "description": "IMAGE-frame/DMSP electron-energy model comparison",
        "source_detector_crossings": str(matches_path.resolve()),
        "source_annotations": str(annotations_path.resolve()),
        "source_kp": str(kp_path.resolve()),
        "annotation_status": "-1 absent, 0 rejected, 1 accepted, 2 deferred",
        "dmsp_energy_reduction": "sum(Q) / sum(Q / E)",
        "dmsp_energy_fractional_std_maximum": MAX_ENERGY_FRACTIONAL_STD,
        "dmsp_flux_fractional_std_maximum": MAX_FLUX_FRACTIONAL_STD,
        "temporal_support_seconds": ", ".join(map(str, support_seconds)),
        "temporal_support_convention": "[image_time - L/2, image_time + L/2)",
        "zp_mlt_interpolation_step_hours": ZP_MLT_STEP_HOURS,
        "zp_note": "grid-independent latitude collapse evaluated by Kp and MLT",
        "hardy_source": "Hardy et al. (1987), Tables 4/5; doi:10.1029/JA092iA11p12275",
        "hardy_energy_reduction": "sum(model energy flux) / sum(model number flux) on valid DMSP samples",
        "hardy_note": "uncalibrated full spatial model; corrected-geomagnetic model evaluated on supplied DMSP coordinates; full-spectrum moments, extrapolated to 100 keV"
    })

    units = {
        "dmsp_energy": "keV",
        "dmsp_flux": "eV cm-2 sr-1 s-1",
        "mlt": "hour",
        "mlat": "degree",
        "detector_separation_median_deg": "degree",
        "detector_separation_maximum_deg": "degree",
        "image_ratio_energy_median": "keV",
        "image_ratio_energy_std_median": "keV",
        "zp_pointwise_energy": "keV",
        "zp_pointwise_flux": "mW m-2",
        "hardy_pointwise_energy": "keV",
        "hardy_pointwise_flux": "keV cm-2 sr-1 s-1",
        "zp_area_mean_energy": "keV",
        "zp_area_median_energy": "keV",
        "zp_flux_weighted_energy": "keV",
        "zp_latitude_spread": "keV",
        "zp_selected_lower_mlat": "degree",
        "zp_selected_upper_mlat": "degree"
    }
    for name, unit in units.items():
        if name in dataset:
            dataset[name].attrs["units"] = unit

    encoding = {
        name: {"zlib": True, "complevel": 4}
        for name in dataset.data_vars
        if dataset[name].dtype.kind not in ("O", "U", "M")
    }
    filename.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_netcdf(filename, engine="netcdf4", encoding=encoding)
    dataset.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matches", type=Path, default=DEFAULT_MATCHES)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--kp", type=Path, default=DEFAULT_KP_PATH)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--support-seconds", type=int, nargs="+",
        default=list(DEFAULT_SUPPORT_SECONDS)
    )
    parser.add_argument("--orbit", type=int, nargs="+")
    args = parser.parse_args()

    if any(seconds <= 0 for seconds in args.support_seconds):
        raise ValueError("support seconds must be positive")
    support_seconds = tuple(sorted(set(args.support_seconds)))

    files = sorted(args.matches.glob("or_*.nc"))
    if args.orbit:
        selected = set(args.orbit)
        files = [
            filename for filename in files
            if int(filename.stem.split("_")[-1]) in selected
        ]
    if not files:
        raise FileNotFoundError(f"No detector crossing files found in {args.matches}")

    annotation_status = load_annotation_status(args.annotations)
    kp_series = load_gfz_kp(args.kp)
    records = []
    for filename in tqdm(files, desc="Reduce detector crossings"):
        records.extend(reduce_orbit(
            filename, support_seconds, annotation_status, kp_series
        ))
    if not records:
        raise ValueError("No DMSP crossing records were produced")

    comparison = pd.DataFrame(records)
    comparison = add_zhang_paxton(comparison, args.kp)
    save_comparison(
        comparison, args.output, args.matches, args.annotations, args.kp,
        support_seconds
    )

    accepted = comparison["annotation_status"].eq(1).sum()
    print(f"Saved {len(comparison):,} records to {args.output}")
    print(f"Accepted annotation records: {accepted:,}")


if __name__ == "__main__":
    main()
