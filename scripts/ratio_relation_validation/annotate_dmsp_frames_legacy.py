"""Accept or reject IMAGE frames that contain matched DMSP measurements."""

#%% Imports
import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from icreader import load as icload

#%% Paths and saved columns
REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_MATCHES = REPOSITORY / "data" / "matches.nc"
DEFAULT_IMAGES = Path("/home/bing/Dropbox/work/data/IMAGE_FUV/precipitation_IR_P2")
DEFAULT_OUTPUT = REPOSITORY / "data" / "dmsp_frame_annotations.csv"

COLUMNS = [
    "orbit", "frame_id", "img_time", "accepted", "satellites",
    "n_dmsp_samples", "annotated_utc",
]

#%% Read and write annotations
def read_completed(filename):
    """Return orbit/timestamp/satellite keys already saved."""

    if not filename.exists():
        return set()

    completed = set()
    with filename.open(newline="") as file:
        for row in csv.DictReader(file):
            satellites = row["satellites"].split(";")
            if len(satellites) == 1:
                completed.add(
                    (int(row["orbit"]), row["img_time"], satellites[0].upper())
                )

    return completed


def save_annotation(filename, annotation):
    """Append one decision immediately so the session can resume later."""

    filename.parent.mkdir(parents=True, exist_ok=True)
    write_header = not filename.exists()
    with filename.open("a", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=COLUMNS)
        if write_header:
            writer.writeheader()
        writer.writerow(annotation)

#%% Display one matched frame
def colour_limit(values):
    finite = values[np.isfinite(values) & (values >= 0)]
    return max(float(np.nanpercentile(finite, 99)), 1.0) if finite.size else 1.0


def plot_frame(image, orbit, frame, matches, active_satellite, progress):
    """Display one active DMSP track and other tracks for context."""

    fields = [
        (np.asarray(image.wic_corrected[frame]), "Corrected WIC"),
        (np.asarray(image.si13_corrected[frame]), "Corrected SI13"),
    ]
    figure, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)

    for axis, (values, title) in zip(axes, fields):
        plotted = axis.pcolormesh(
            image.grid.xi, image.grid.eta, values,
            shading="auto", cmap="viridis", vmin=0,
            vmax=colour_limit(values),
        )

        for satellite in np.unique(matches.dmsp_sat.values):
            track = matches.where(matches.dmsp_sat == satellite, drop=True)
            track = track.sortby("dmsp_time")
            satellite_name = str(satellite).upper()
            xi, eta = image.grid.projection.geo2cube(
                track.dmsp_mlt.values * 15,
                track.dmsp_mlat.values,
                set_points_off_cube_to_nan=True,
            )

            if satellite_name == active_satellite:
                axis.plot(xi, eta, color="black", linewidth=7)
                axis.plot(xi, eta, color="magenta", linewidth=4,
                          label=f"{satellite_name} — annotate")
                axis.plot(xi[0], eta[0], "o", color="magenta",
                          markeredgecolor="black", markersize=9)
                axis.plot(xi[-1], eta[-1], "^", color="magenta",
                          markeredgecolor="black", markersize=10)
            else:
                axis.plot(xi, eta, color="0.55", linewidth=2,
                          label=f"{satellite_name} — other")

        axis.set_title(title)
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        axis.legend(loc="lower left")
        figure.colorbar(plotted, ax=axis, shrink=0.8, label="counts")

    image_time = np.datetime_as_string(
        np.asarray(image.time[frame], dtype="datetime64[ns]"), unit="s"
    )
    active_matches = matches.where(
        matches.dmsp_sat.astype(str).str.upper() == active_satellite,
        drop=True,
    )
    figure.suptitle(
        f"Orbit {orbit:04d}, frame {frame:03d}, {image_time} UTC\n"
        f"Annotating {active_satellite}: {active_matches.sizes['index']} samples"
        f" — {progress}\n"
        "Press 0 to reject, 1 to accept, or q to quit"
    )
    return figure, image_time


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

def annotate(matches_file, image_dir, output, selected_orbits=None):
    """Step through each unique orbit/frame represented in matches.nc."""

    completed = read_completed(output)

    with xr.open_dataset(matches_file) as matches:
        # Loading only the orbit column avoids putting the full 2.8-GB file in memory.
        orbit_values = np.asarray(matches.orbit.values)
        orbits, starts, counts = np.unique(
            orbit_values, return_index=True, return_counts=True
        )
        if selected_orbits is not None:
            selected_orbits = set(selected_orbits)

        ranges = [
            (int(orbit), int(start), int(start + count))
            for orbit, start, count in zip(orbits, starts, counts)
            if selected_orbits is None or int(orbit) in selected_orbits
        ]

        for orbit, start, stop in ranges:
            image_file = image_dir / f"or_{orbit:04d}.nc"
            if not image_file.exists():
                print(f"Missing image product: {image_file}")
                continue

            orbit_matches = matches.isel(index=slice(start, stop)).load()
            frames = np.unique(orbit_matches.frame_id.values).astype(int)
            image = None

            for number, frame in enumerate(frames, start=1):
                frame_matches = orbit_matches.where(
                    orbit_matches.frame_id == frame, drop=True
                )
                match_time = np.datetime_as_string(
                    frame_matches.img_time.values[0], unit="s"
                )
                satellite_values = sorted(
                    np.unique(frame_matches.dmsp_sat.values), key=str
                )

                for satellite_number, satellite_value in enumerate(
                    satellite_values, start=1
                ):
                    satellite = str(satellite_value).upper()
                    if (orbit, match_time, satellite) in completed:
                        continue

                    satellite_matches = frame_matches.where(
                        frame_matches.dmsp_sat == satellite_value, drop=True
                    )

                    # A ratio is usable only where WIC and SI13 are both finite.
                    wic = np.asarray(satellite_matches.img_wic.values)
                    si13 = np.asarray(satellite_matches.img_s13.values)
                    has_image_data = np.any(
                        np.isfinite(wic) & np.isfinite(si13)
                    )

                    dmsp_mlt = np.asarray(satellite_matches.dmsp_mlt.values)
                    has_nightside_data = np.any(
                        np.isfinite(dmsp_mlt)
                        & ((dmsp_mlt >= 18) | (dmsp_mlt <= 6))
                    )

                    if not has_image_data:
                        decision = "0"
                        image_time = match_time
                        print(
                            f"Auto-rejected orbit {orbit:04d}, "
                            f"frame {frame:03d}, {satellite}: "
                            "no matched WIC/SI13 data"
                        )
                    elif not has_nightside_data:
                        decision = "2"
                        image_time = match_time
                        print(
                            f"Deferred orbit {orbit:04d}, frame {frame:03d}, "
                            f"{satellite}: DMSP samples only the dayside"
                        )
                    else:
                        if image is None:
                            image = icload(image_file)

                        progress = (
                            f"frame {number}/{len(frames)}, satellite "
                            f"{satellite_number}/{len(satellite_values)}"
                        )
                        figure, image_time = plot_frame(
                            image, orbit, frame, frame_matches,
                            satellite, progress,
                        )
                        if image_time != match_time:
                            plt.close(figure)
                            raise ValueError(
                                f"orbit {orbit:04d} frame {frame} has time "
                                f"{image_time}, but matches.nc has {match_time}"
                            )

                        decision = read_keypress(figure)
                        if decision == "q":
                            return

                    save_annotation(output, {
                        "orbit": orbit,
                        "frame_id": frame,
                        "img_time": image_time,
                        "accepted": int(decision),
                        "satellites": satellite,
                        "n_dmsp_samples": satellite_matches.sizes["index"],
                        "annotated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    })
                    completed.add((orbit, image_time, satellite))


#%% Command line

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matches", type=Path, default=DEFAULT_MATCHES)
    parser.add_argument("--image-dir", type=Path, default=DEFAULT_IMAGES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--orbit", type=int, nargs="+",
                        help="Annotate only these orbit numbers")
    args = parser.parse_args()
    annotate(args.matches, args.image_dir, args.output, args.orbit)


if __name__ == "__main__":
    main()
