"""Generate lightweight contact sheets for manual FUV frame-quality review.

The script reads the original IDL frames with the same quality-relevant loader
settings used by the production background pipeline, but stops before
``backgroundmodel_BS`` removes quality-0/1 frames. Each sensor/orbit sheet
contains the unique first and last 18 northern frames in a 6-by-6 layout. A
per-sheet CSV and a combined manifest retain the stable source identity and
panel location needed by the local annotation script.
"""

#%% Imports

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import fuvpy as fuv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#%% Configuration

REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_BASE = REPOSITORY / "example_data"
PANELS_PER_END = 18
PANEL_ROWS = 6
PANEL_COLUMNS = 6

SENSOR_SETTINGS = {
    "wic": {
        "label": "WIC", "index": "wicfiles.h5",
        "folder": "wic_data", "reflat": True
    },
    "si12": {
        "label": "SI12", "index": "s12files.h5",
        "folder": "s12_data", "reflat": False
    },
    "si13": {
        "label": "SI13", "index": "s13files.h5",
        "folder": "s13_data", "reflat": False
    }
}
QUALITY_COLOURS = {0: "red", 1: "#e6b800", 2: "#00a651"}
MANIFEST_COLUMNS = [
    "sheet", "panel_index", "panel_row", "panel_column",
    "panel_left", "panel_top", "panel_right", "panel_bottom",
    "orbit", "sensor", "source_file", "source_index",
    "northern_frame_index", "time", "automatic_quality",
    "mapped_fraction", "maximum_absolute_mlat",
    "median_row_relative_spread", "turn_on_reference_correlation",
    "hv_mcp", "hv_phos"
]


#%% Selection and file helpers

def select_edge_frames(frame_count, panels_per_end=PANELS_PER_END):
    """Return unique chronological indices from both ends of one orbit."""

    if frame_count < 0:
        raise ValueError("frame_count must not be negative")
    if panels_per_end < 1:
        raise ValueError("panels_per_end must be at least one")

    first = np.arange(min(panels_per_end, frame_count), dtype=int)
    last_start = max(frame_count - panels_per_end, 0)
    last = np.arange(last_start, frame_count, dtype=int)
    return np.unique(np.concatenate([first, last]))


def atomic_write_csv(table, filename, columns=None):
    """Replace one CSV only after its complete temporary file is written."""

    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    temporary = filename.with_name(filename.name + ".partial")
    table.to_csv(temporary, index=False, columns=columns)
    temporary.replace(filename)


def read_sensor_index(base, sensor):
    """Read and normalize one production HDF5 orbit index."""

    filename = base / SENSOR_SETTINGS[sensor]["index"]
    files = pd.read_hdf(filename, key="data").reset_index()
    required = {"filename", "orbit"}
    missing = required.difference(files.columns)
    if missing:
        raise ValueError(f"{filename} is missing {sorted(missing)}")
    files["orbit"] = files["orbit"].astype(int)
    return files


def source_files_for_orbit(base, sensor, files, orbit):
    """Return lexically sorted raw files for one indexed sensor orbit."""

    folder = base / SENSOR_SETTINGS[sensor]["folder"]
    names = files.loc[files["orbit"] == orbit, "filename"].astype(str)
    paths = sorted(folder / name for name in names)
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing raw input: {missing[0]}")
    return paths


#%% Raw frame loading

def read_orbit_frames(source_files, sensor):
    """Load pre-background frames and retain the production quality flags."""

    filenames = [str(path) for path in source_files]
    # Measurement variance and reconstructed viewing geometry are downstream
    # fields. They do not enter fuvpy's frame-quality calculation, so omitting
    # them makes this review pass substantially cheaper without changing flags.
    images = fuv.read_idl(
        filenames,
        dzalim=75,
        reflat=SENSOR_SETTINGS[sensor]["reflat"],
        measurement_variance=False,
        viewing_geometry=False
    )
    if images.sizes.get("date", 0) != len(source_files):
        raise ValueError("read_idl returned a different number of frames")

    images = images.assign_coords(
        source_index=("date", np.arange(len(source_files), dtype=int)),
        source_file=("date", [path.name for path in source_files])
    )
    north = np.asarray(images.hemisphere.values).astype(str) == "north"
    images = images.isel(date=np.flatnonzero(north))
    images = images.assign_coords(
        northern_frame_index=("date", np.arange(images.sizes["date"], dtype=int))
    )
    return images


def frame_metrics(images):
    """Calculate compact diagnostics already used by the quality classifier."""

    image = np.asarray(images.img.values, dtype=float)
    mlat = np.asarray(images.mlat.values, dtype=float)
    sza = np.asarray(images.sza.values, dtype=float)
    dza = np.asarray(images.dza.values, dtype=float)
    mapped = (
        np.isfinite(image) & np.isfinite(mlat)
        & np.isfinite(sza) & np.isfinite(dza)
    )
    mapped_fraction = np.mean(mapped, axis=(1, 2))

    maximum_absolute_mlat = np.full(images.sizes["date"], np.nan)
    for frame in range(images.sizes["date"]):
        if np.any(mapped[frame]):
            maximum_absolute_mlat[frame] = np.max(
                np.abs(mlat[frame][mapped[frame]])
            )

    return {
        "mapped_fraction": mapped_fraction,
        "maximum_absolute_mlat": maximum_absolute_mlat,
        "median_row_relative_spread": np.asarray(
            images.median_row_relative_spread.values, dtype=float
        ),
        "turn_on_reference_correlation": np.asarray(
            images.turn_on_reference_correlation.values, dtype=float
        ),
        "hv_mcp": np.asarray(images.hv_mcp.values, dtype=float),
        "hv_phos": np.asarray(images.hv_phos.values, dtype=float)
    }


#%% Contact-sheet rendering

def common_colour_limit(values):
    """Return one robust nonnegative scale shared by every selected panel."""

    finite = values[np.isfinite(values) & (values >= 0)]
    if not finite.size:
        return 1.0
    return max(float(np.nanpercentile(finite, 99.5)), 1.0)


def panel_bounds(axis):
    """Return one subplot rectangle in top-left normalized image coordinates."""

    left, bottom, width, height = axis.get_position().bounds
    return left, 1 - bottom - height, left + width, 1 - bottom


def make_contact_sheet(images, orbit, sensor, output_file, dpi=150):
    """Write one 6-by-6 contact sheet and return its panel manifest."""

    frame_count = images.sizes.get("date", 0)
    selected = select_edge_frames(frame_count)
    if not selected.size:
        return pd.DataFrame(columns=MANIFEST_COLUMNS)

    selected_images = np.asarray(images.img.isel(date=selected).values, dtype=float)
    colour_limit = common_colour_limit(selected_images)
    qualities = np.asarray(images.frame_quality.values, dtype=int)
    metrics = frame_metrics(images)

    figure, axes = plt.subplots(
        PANEL_ROWS, PANEL_COLUMNS, figsize=(13, 13), constrained_layout=False
    )
    figure.subplots_adjust(
        left=0.025, right=0.99, bottom=0.025, top=0.94,
        wspace=0.08, hspace=0.20
    )
    flat_axes = axes.ravel()

    rows = []
    relative_sheet = Path(sensor) / output_file.name
    for panel_index, frame_index in enumerate(selected):
        axis = flat_axes[panel_index]
        quality = int(qualities[frame_index])
        axis.imshow(
            selected_images[panel_index], origin="lower", cmap="magma",
            vmin=0, vmax=colour_limit, interpolation="nearest",
            aspect="equal"
        )
        axis.set_xticks([])
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_visible(True)
            spine.set_color(QUALITY_COLOURS[quality])
            spine.set_linewidth(4)

        timestamp = np.datetime64(images.date.values[frame_index], "ns")
        time_text = np.datetime_as_string(timestamp, unit="s")
        axis.set_title(
            f"{int(images.northern_frame_index.values[frame_index]):03d}  "
            f"{time_text[11:19]}  Q{quality}",
            fontsize=6, pad=2
        )

        left, top, right, bottom = panel_bounds(axis)
        row = {
            "sheet": relative_sheet.as_posix(),
            "panel_index": panel_index,
            "panel_row": panel_index // PANEL_COLUMNS,
            "panel_column": panel_index % PANEL_COLUMNS,
            "panel_left": left,
            "panel_top": top,
            "panel_right": right,
            "panel_bottom": bottom,
            "orbit": int(orbit),
            "sensor": SENSOR_SETTINGS[sensor]["label"],
            "source_file": str(images.source_file.values[frame_index]),
            "source_index": int(images.source_index.values[frame_index]),
            "northern_frame_index": int(
                images.northern_frame_index.values[frame_index]
            ),
            "time": time_text,
            "automatic_quality": quality
        }
        for name, values in metrics.items():
            row[name] = float(values[frame_index])
        rows.append(row)

    for axis in flat_axes[selected.size:]:
        axis.set_visible(False)

    counts = np.bincount(qualities, minlength=3)
    figure.suptitle(
        f"{SENSOR_SETTINGS[sensor]['label']} orbit {orbit:04d}: "
        f"first/last {PANELS_PER_END} northern frames\n"
        f"all northern frames={frame_count}; quality 0/1/2="
        f"{counts[0]}/{counts[1]}/{counts[2]}; display 0--{colour_limit:.3g} counts",
        fontsize=12
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_file, dpi=dpi)
    plt.close(figure)
    return pd.DataFrame(rows, columns=MANIFEST_COLUMNS)


#%% One sensor/orbit unit and restart handling

def process_orbit(base, output, sensor, orbit, files, dpi=150, overwrite=False):
    """Generate one restartable PNG/CSV review unit."""

    sensor_output = output / sensor
    image_file = sensor_output / f"or_{orbit:04d}.png"
    manifest_file = sensor_output / f"or_{orbit:04d}.csv"
    if image_file.exists() and manifest_file.exists() and not overwrite:
        return "skipped", manifest_file

    source_files = source_files_for_orbit(base, sensor, files, orbit)
    if not source_files:
        return "empty", None
    images = read_orbit_frames(source_files, sensor)
    if images.sizes.get("date", 0) == 0:
        return "empty", None

    manifest = make_contact_sheet(
        images, orbit, sensor, image_file, dpi=dpi
    )
    atomic_write_csv(manifest, manifest_file, MANIFEST_COLUMNS)
    return "written", manifest_file


def rebuild_combined_manifest(output):
    """Combine completed per-sheet manifests after an interrupted-safe run."""

    files = sorted(
        path for sensor in SENSOR_SETTINGS
        for path in (output / sensor).glob("or_*.csv")
    )
    if files:
        manifest = pd.concat(
            [pd.read_csv(filename) for filename in files], ignore_index=True
        )
        manifest = manifest.sort_values(
            ["sensor", "orbit", "panel_index"]
        ).reset_index(drop=True)
    else:
        manifest = pd.DataFrame(columns=MANIFEST_COLUMNS)
    atomic_write_csv(manifest, output / "manifest.csv", MANIFEST_COLUMNS)
    return manifest


#%% Command line

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-input", type=Path, default=DEFAULT_BASE,
        help="directory containing the sensor HDF5 indices and *_data folders"
    )
    parser.add_argument(
        "--output", type=Path,
        help="review-package directory; defaults to BASE_INPUT/frame_quality_review"
    )
    parser.add_argument(
        "--sensor", nargs="+", choices=tuple(SENSOR_SETTINGS),
        default=list(SENSOR_SETTINGS), help="sensors to render"
    )
    parser.add_argument("--orbit", type=int, nargs="+")
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    base = args.base_input.expanduser().resolve()
    output = (
        args.output.expanduser().resolve()
        if args.output is not None else base / "frame_quality_review"
    )
    output.mkdir(parents=True, exist_ok=True)

    summary = {"written": 0, "skipped": 0, "empty": 0, "failed": 0}
    failures = []
    selected_orbits = None if args.orbit is None else set(args.orbit)

    for sensor in args.sensor:
        files = read_sensor_index(base, sensor)
        orbits = sorted(files.orbit.unique())
        if selected_orbits is not None:
            orbits = [orbit for orbit in orbits if orbit in selected_orbits]

        for position, orbit in enumerate(orbits, start=1):
            print(
                f"{SENSOR_SETTINGS[sensor]['label']} orbit {orbit:04d} "
                f"({position}/{len(orbits)})",
                flush=True
            )
            try:
                status, _ = process_orbit(
                    base, output, sensor, int(orbit), files,
                    dpi=args.dpi, overwrite=args.overwrite
                )
            except Exception as error:
                status = "failed"
                failures.append({
                    "sensor": SENSOR_SETTINGS[sensor]["label"],
                    "orbit": int(orbit), "error": str(error)
                })
                print(f"  failed: {error}", flush=True)
            summary[status] += 1

    manifest = rebuild_combined_manifest(output)
    atomic_write_csv(
        pd.DataFrame(failures, columns=["sensor", "orbit", "error"]),
        output / "failures.csv", ["sensor", "orbit", "error"]
    )
    print(f"Manifest rows: {len(manifest):,}")
    print("Summary:")
    for name, count in summary.items():
        print(f" {name}: {count:,}")


if __name__ == "__main__":
    main()
