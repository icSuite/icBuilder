"""Annotate detector-matched DMSP segments using compact CS images."""

#%% Imports

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path

from icreader import open_product
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


#%% Paths and saved columns

REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_MATCHES = Path(
    "/home/bing/Dropbox/work/temp_storage/icBuilder_pipeline_test/"
    "ratio_relation_validation/detector_crossings"
)
DEFAULT_IMAGES = (
    Path.home() / "dynamit_server/IMAGE_FUV/precipitation_cs/"
    "image_apex_130km_46x46_v1/IR_hardy"
)
DEFAULT_OUTPUT = REPOSITORY / "data/dmsp_frame_annotations.csv"

COLUMNS = [
    "orbit", "frame_id", "img_time", "accepted", "satellites",
    "n_dmsp_samples", "annotated_utc",
]
REQUIRED_MATCH_FIELDS = {
    "image_frame", "image_time", "dmsp_sat", "dmsp_time",
    "dmsp_mlat", "dmsp_mlt", "detector_row", "detector_column",
    "img_wic", "img_si13", "wic_sza",
}


#%% Read and write annotations

def annotation_key(row):
    """Return the stable key used by one single-satellite annotation."""

    satellites = row["satellites"].split(";")
    if len(satellites) != 1:
        return None
    return int(row["orbit"]), row["img_time"], satellites[0].upper()


def read_annotation_statuses(filename):
    """Return the saved status for each single-satellite annotation key."""

    if not filename.exists():
        return {}

    statuses = {}
    with filename.open(newline="") as file:
        for row in csv.DictReader(file):
            key = annotation_key(row)
            if key is not None:
                statuses[key] = int(row["accepted"])
    return statuses


def completed_keys(statuses, revisit_deferred=False):
    """Return saved keys that should remain completed for this run."""

    return {
        key for key, status in statuses.items()
        if not (revisit_deferred and status == 2)
    }


def read_completed(filename, revisit_deferred=False):
    """Return saved keys, optionally leaving status-2 entries unfinished."""

    statuses = read_annotation_statuses(filename)
    return completed_keys(statuses, revisit_deferred=revisit_deferred)


def save_annotation(filename, annotation, replace_deferred=False):
    """Save one decision, replacing its old status-2 row when requested."""

    filename.parent.mkdir(parents=True, exist_ok=True)

    if replace_deferred and filename.exists():
        with filename.open(newline="") as file:
            rows = list(csv.DictReader(file))

        key = annotation_key(annotation)
        for index, row in enumerate(rows):
            if annotation_key(row) == key and int(row["accepted"]) == 2:
                rows[index] = annotation
                temporary = filename.with_name(filename.name + ".partial")
                with temporary.open("w", newline="") as file:
                    writer = csv.DictWriter(file, fieldnames=COLUMNS)
                    writer.writeheader()
                    writer.writerows(rows)
                temporary.replace(filename)
                return True

    write_header = not filename.exists()
    with filename.open("a", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=COLUMNS)
        if write_header:
            writer.writeheader()
        writer.writerow(annotation)
    return False


#%% Candidate selection

def finite_values(values):
    values = np.asarray(values, dtype=float)
    return values[np.isfinite(values)]


def has_minimum_sza(matches, minimum_sza):
    """Return whether any matched sample reaches the requested darkness."""

    sza = finite_values(matches.wic_sza.values)
    return bool(sza.size and np.any(sza >= minimum_sza))


def has_image_data(matches):
    """Return whether WIC and SI13 are simultaneously finite on the track."""

    wic = np.asarray(matches.img_wic.values, dtype=float)
    si13 = np.asarray(matches.img_si13.values, dtype=float)
    return bool(np.any(np.isfinite(wic) & np.isfinite(si13)))


def has_nightside_mlt(matches):
    """Preserve the legacy 18--06 MLT proxy when no SZA filter is requested."""

    mlt = np.asarray(matches.dmsp_mlt.values, dtype=float)
    return bool(np.any(np.isfinite(mlt) & ((mlt >= 18) | (mlt <= 6))))


def orbit_number(filename):
    return int(filename.stem.split("_")[-1])


def find_match_files(matches_dir, selected_orbits=None):
    """Find per-orbit detector crossing files in numerical order."""

    files = sorted(matches_dir.glob("or_*.nc"), key=orbit_number)
    if selected_orbits is not None:
        selected = set(selected_orbits)
        files = [filename for filename in files if orbit_number(filename) in selected]
    if not files:
        raise FileNotFoundError(f"No detector crossing files found in {matches_dir}")
    return files


def validate_crossing_file(matches, filename):
    """Check the small part of the crossing contract used by annotation."""

    if matches.sizes.get("sample", 0) == 0:
        return
    if matches.attrs.get("product_type") != "dmsp_image_crossings":
        raise ValueError(f"{filename} is not a DMSP crossing product")
    if matches.attrs.get("representation") != "detector":
        raise ValueError(f"{filename} is not detector crossing data")
    missing = REQUIRED_MATCH_FIELDS.difference(matches.variables)
    if missing:
        raise ValueError(f"{filename} is missing {sorted(missing)}")


#%% Display one matched frame

def colour_limit(values):
    finite = values[np.isfinite(values) & (values >= 0)]
    return max(float(np.nanpercentile(finite, 99)), 1.0) if finite.size else 1.0


def cs_frame_for_time(image, match_time, suggested_frame):
    """Resolve a crossing frame against the CS product by absolute time."""

    times = np.asarray(image.time, dtype="datetime64[ns]")
    target = np.datetime64(match_time, "ns")
    if (
        0 <= suggested_frame < times.size
        and times[suggested_frame] == target
    ):
        return suggested_frame

    found = np.flatnonzero(times == target)
    if found.size != 1:
        raise ValueError(
            f"CS product contains {found.size} frames at "
            f"{np.datetime_as_string(target, unit='s')}"
        )
    return int(found[0])


def plot_frame(
    image, orbit, crossing_frame, cs_frame, match_time, matches,
    active_satellite, progress,
):
    """Display one active detector crossing over compact CS WIC and SI13."""

    fields = [
        (np.asarray(image.wic_corrected[cs_frame]), "Corrected WIC"),
        (np.asarray(image.si13_corrected[cs_frame]), "Corrected SI13"),
    ]
    figure, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)

    satellite_values = np.asarray(matches.dmsp_sat.values).astype(str)
    for axis, (values, title) in zip(axes, fields):
        plotted = axis.pcolormesh(
            image.grid.xi, image.grid.eta, values,
            shading="auto", cmap="viridis", vmin=0,
            vmax=colour_limit(values),
        )

        for satellite_value in np.unique(satellite_values):
            track = matches.where(matches.dmsp_sat == satellite_value, drop=True)
            track = track.sortby("dmsp_time")
            satellite_name = str(satellite_value).upper()
            xi, eta = image.grid.projection.geo2cube(
                track.dmsp_mlt.values * 15,
                track.dmsp_mlat.values,
                set_points_off_cube_to_nan=True,
            )

            if satellite_name == active_satellite:
                axis.plot(xi, eta, color="black", linewidth=7)
                axis.plot(
                    xi, eta, color="magenta", linewidth=4,
                    label=f"{satellite_name} — annotate",
                )
                axis.plot(
                    xi[0], eta[0], "o", color="magenta",
                    markeredgecolor="black", markersize=9,
                )
                axis.plot(
                    xi[-1], eta[-1], "^", color="magenta",
                    markeredgecolor="black", markersize=10,
                )
            else:
                axis.plot(
                    xi, eta, color="0.55", linewidth=2,
                    label=f"{satellite_name} — other",
                )

        axis.set_title(title)
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        axis.legend(loc="lower left")
        figure.colorbar(plotted, ax=axis, shrink=0.8, label="counts")

    active = np.char.upper(satellite_values) == active_satellite
    active_sza = finite_values(matches.wic_sza.values[active])
    sza_text = (
        f"WIC SZA {active_sza.min():.1f}--{active_sza.max():.1f}°"
        if active_sza.size else "WIC SZA unavailable"
    )
    image_time = np.datetime_as_string(
        np.datetime64(match_time, "ns"), unit="s"
    )
    figure.suptitle(
        f"Orbit {orbit:04d}, detector frame {crossing_frame:03d}, "
        f"{image_time} UTC\n"
        f"Annotating {active_satellite}: {active.sum()} samples, {sza_text}"
        f" — {progress}\n"
        "Press 0 to reject, 1 to accept, or q to quit"
    )
    return figure


def read_keypress(figure):
    """Wait for one annotation key from the active figure window."""

    decision = []

    def key_pressed(event):
        if event.key in ("0", "1", "q"):
            decision.append(event.key)
            plt.close(figure)

    figure.canvas.mpl_connect("key_press_event", key_pressed)
    plt.show()
    return decision[0] if decision else "q"


#%% Annotation loop

def open_cs_image(filename):
    """Open and identify the compact visual-context product."""

    image = open_product(filename)
    if image.product_type != "precipitation_cs" or image.representation != "cs":
        image.close()
        raise ValueError(f"{filename} is not CS Product 2")
    return image


def annotate(
    matches_dir, image_dir, output, selected_orbits=None, minimum_sza=None,
    revisit_deferred=False,
):
    """Step through orbit/frame/satellite segments not already annotated."""

    if minimum_sza is not None and not np.isfinite(minimum_sza):
        raise ValueError("minimum SZA must be finite")

    statuses = read_annotation_statuses(output)
    completed = completed_keys(statuses, revisit_deferred=revisit_deferred)
    files = find_match_files(matches_dir, selected_orbits)
    summary = {
        "completed": 0,
        "non_dark_skipped": 0,
        "unusable_deferred_skipped": 0,
        "detector_pair_missing_skipped": 0,
        "mlt_deferred": 0,
        "deferred_replaced": 0,
        "presented": 0,
    }

    for matches_file in files:
        orbit = orbit_number(matches_file)
        with xr.open_dataset(matches_file) as opened:
            validate_crossing_file(opened, matches_file)
            if opened.sizes.get("sample", 0) == 0:
                continue
            orbit_matches = opened.load()

        frames = np.unique(orbit_matches.image_frame.values).astype(int)
        image = None
        try:
            for frame_number, frame in enumerate(frames, start=1):
                frame_matches = orbit_matches.where(
                    orbit_matches.image_frame == frame, drop=True
                )
                frame_times = np.unique(
                    np.asarray(frame_matches.image_time.values, dtype="datetime64[ns]")
                )
                if frame_times.size != 1:
                    raise ValueError(
                        f"orbit {orbit:04d} frame {frame} has "
                        f"{frame_times.size} IMAGE times"
                    )
                match_time = frame_times[0]
                match_time_text = np.datetime_as_string(match_time, unit="s")
                satellite_values = sorted(
                    np.unique(frame_matches.dmsp_sat.values), key=str
                )

                for satellite_number, satellite_value in enumerate(
                    satellite_values, start=1
                ):
                    satellite = str(satellite_value).upper()
                    key = (orbit, match_time_text, satellite)
                    if key in completed:
                        summary["completed"] += 1
                        continue

                    satellite_matches = frame_matches.where(
                        frame_matches.dmsp_sat == satellite_value, drop=True
                    )

                    # The explicit SZA option is a transient filter. Nothing is
                    # saved, so this candidate remains available in a later run.
                    if (
                        minimum_sza is not None
                        and not has_minimum_sza(satellite_matches, minimum_sza)
                    ):
                        summary["non_dark_skipped"] += 1
                        continue

                    if (
                        revisit_deferred
                        and statuses.get(key) == 2
                        and not has_image_data(satellite_matches)
                    ):
                        summary["unusable_deferred_skipped"] += 1
                        continue
                    elif not has_image_data(satellite_matches):
                        summary["detector_pair_missing_skipped"] += 1
                        print(
                            f"Skipped orbit {orbit:04d}, frame {frame:03d}, "
                            f"{satellite}: no jointly finite detector-matched "
                            "WIC/SI13 samples; annotation unchanged"
                        )
                        continue
                    elif minimum_sza is None and not has_nightside_mlt(
                        satellite_matches
                    ):
                        decision = "2"
                        summary["mlt_deferred"] += 1
                        print(
                            f"Deferred orbit {orbit:04d}, frame {frame:03d}, "
                            f"{satellite}: DMSP samples only the dayside"
                        )
                    else:
                        image_file = image_dir / f"or_{orbit:04d}.nc"
                        if not image_file.exists():
                            raise FileNotFoundError(
                                f"Missing CS image product: {image_file}"
                            )
                        if image is None:
                            image = open_cs_image(image_file)
                        cs_frame = cs_frame_for_time(image, match_time, frame)
                        progress = (
                            f"frame {frame_number}/{len(frames)}, satellite "
                            f"{satellite_number}/{len(satellite_values)}"
                        )
                        figure = plot_frame(
                            image, orbit, frame, cs_frame, match_time,
                            frame_matches, satellite, progress,
                        )
                        decision = read_keypress(figure)
                        if decision == "q":
                            return summary
                        summary["presented"] += 1

                    replaced = save_annotation(output, {
                        "orbit": orbit,
                        "frame_id": frame,
                        "img_time": match_time_text,
                        "accepted": int(decision),
                        "satellites": satellite,
                        "n_dmsp_samples": satellite_matches.sizes["sample"],
                        "annotated_utc": datetime.now(timezone.utc).isoformat(
                            timespec="seconds"
                        ),
                    }, replace_deferred=revisit_deferred)
                    if replaced:
                        summary["deferred_replaced"] += 1
                    statuses[key] = int(decision)
                    completed.add(key)
        finally:
            if image is not None:
                image.close()

    return summary


#%% Command line

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matches", type=Path, default=DEFAULT_MATCHES,
        help="directory containing per-orbit detector crossing files",
    )
    parser.add_argument(
        "--image-dir", type=Path, default=DEFAULT_IMAGES,
        help="directory containing compact precipitation_cs orbit files",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--orbit", type=int, nargs="+", help="annotate only these orbit numbers"
    )
    parser.add_argument(
        "--minimum-sza", type=float,
        help=(
            "show only segments containing at least one sample at or above "
            "this WIC SZA; skipped segments are not written to the CSV"
        ),
    )
    parser.add_argument(
        "--revisit-deferred", action="store_true",
        help=(
            "show existing accepted=2 segments again while keeping accepted=0/1 "
            "completed; each new decision replaces the old status-2 row"
        ),
    )
    args = parser.parse_args()

    summary = annotate(
        args.matches,
        args.image_dir,
        args.output,
        selected_orbits=args.orbit,
        minimum_sza=args.minimum_sza,
        revisit_deferred=args.revisit_deferred,
    )
    print("Annotation summary:")
    for name, count in summary.items():
        print(f" {name}: {count:,}")


if __name__ == "__main__":
    main()
