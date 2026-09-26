"""Compare Product-2 smoothing techniques frame by frame."""

#%% Imports

import argparse
from contextlib import ExitStack
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from icreader import open_product
from tqdm import tqdm


REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_BASE = Path(
    "/home/bing/Dropbox/work/temp_storage/icBuilder_pipeline_test"
)
DEFAULT_OUTPUT = REPOSITORY / "figures/debugging/precipitation_smoothing"
RATIO_MAXIMUM = 136.49


#%% Product discovery and labels

def orbit_number(filename):
    """Return the integer suffix from an or_XXXX.nc filename."""

    text = filename.stem.removeprefix("or_")
    return int(text) if text.isdigit() else None


def discover_product_directories(base, supplied):
    """Use explicit Product-2 directories or discover Hardy test variants."""

    if supplied:
        directories = [path.expanduser() for path in supplied]
    else:
        root = base / "precipitation_detector"
        directories = [root / "IR_hardy"]
        directories.extend(sorted(root.glob("IR_hardy_smooth_*")))

    missing = [directory for directory in directories if not directory.is_dir()]
    if missing:
        raise FileNotFoundError(f"Product-2 directories not found: {missing}")
    return directories


def common_orbits(directories, selected):
    """Find orbit files present in every requested technique directory."""

    available = []
    for directory in directories:
        orbits = {
            orbit_number(filename)
            for filename in directory.glob("or_*.nc")
        }
        available.append({orbit for orbit in orbits if orbit is not None})

    common = set.intersection(*available)
    if selected:
        missing = set(selected).difference(common)
        if missing:
            raise FileNotFoundError(
                f"orbits missing from one or more techniques: {sorted(missing)}"
            )
        return sorted(set(selected))
    return sorted(common)


def technique_label(product):
    """Make a compact row label from Product-2 smoothing provenance."""

    kernel = product.attrs.get("spatial_smoothing_kernel", "none")
    if kernel == "none":
        return "No smoothing"

    wic_width = float(product.attrs["wic_smoothing_width_pixels"])
    si13_width = float(product.attrs["si13_smoothing_width_pixels"])
    if kernel == "gaussian":
        return f"Gaussian: WIC σ={wic_width:g}, SI13 σ={si13_width:g}"
    return f"Boxcar: WIC {wic_width:g}×{wic_width:g}, SI13 {si13_width:g}×{si13_width:g}"


#%% Frame plotting

def positive_limit(images, percentile=99.5):
    """Return one robust positive colour maximum for a sensor column."""

    pieces = [
        image[np.isfinite(image)] for image in images
        if np.any(np.isfinite(image))
    ]
    if not pieces:
        return 1.0
    finite = np.concatenate(pieces)
    positive = finite[finite > 0]
    return float(np.percentile(positive, percentile)) if positive.size else 1.0


def plot_frame(orbit, frame, time, techniques, fuv, output, ratio_maximum):
    """Plot one detector frame with one row per smoothing technique."""

    original = {
        "wic": np.asarray(fuv.read("wic_counts", frame), dtype=float),
        "si13": np.asarray(fuv.read("si13_counts", frame), dtype=float),
    }
    rows = []
    for label, product in techniques:
        if product.spatial_smoothing_kernel == "none":
            wic = original["wic"]
            si13 = original["si13"]
        else:
            wic = np.asarray(product.read("wic_smoothed", frame), dtype=float)
            si13 = np.asarray(product.read("si13_smoothed", frame), dtype=float)
        si12 = np.asarray(product.read("si12", frame), dtype=float)
        ratio = np.asarray(product.read("R", frame), dtype=float)
        rows.append((label, wic, si12, si13, ratio))

    limits = [positive_limit([row[column] for row in rows]) for column in (1, 2, 3)]
    figure, axes = plt.subplots(
        len(rows), 4,
        figsize=(14, max(3.1 * len(rows), 4)),
        squeeze=False,
        constrained_layout=True,
    )
    titles = ("WIC dgimg", "SI12 dgimg", "SI13 dgimg", "WIC/SI13 ratio")
    image_handles = [None] * 4

    for row_index, row in enumerate(rows):
        label, *images = row
        for column, image in enumerate(images):
            maximum = ratio_maximum if column == 3 else limits[column]
            image_handles[column] = axes[row_index, column].imshow(
                image,
                origin="lower",
                interpolation="nearest",
                cmap="viridis",
                vmin=0,
                vmax=maximum,
            )
            axes[row_index, column].set_xticks([])
            axes[row_index, column].set_yticks([])
            if row_index == 0:
                axes[row_index, column].set_title(titles[column])
        axes[row_index, 0].set_ylabel(label)

    for column, handle in enumerate(image_handles):
        figure.colorbar(handle, ax=axes[:, column], shrink=0.82, pad=0.01)

    timestamp = np.datetime_as_string(np.datetime64(time), unit="s")
    figure.suptitle(f"Orbit {orbit:04d}, frame {frame:03d}, {timestamp} UTC")
    orbit_output = output / f"or_{orbit:04d}"
    orbit_output.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        orbit_output / f"frame_{frame:03d}.png",
        dpi=160,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_orbit(orbit, directories, output, start_frame, stop_frame, ratio_maximum):
    """Open all techniques and generate every selected frame."""

    with ExitStack() as stack:
        products = [
            stack.enter_context(open_product(directory / f"or_{orbit:04d}.nc"))
            for directory in directories
        ]
        techniques = sorted(
            [(technique_label(product), product) for product in products],
            key=lambda item: item[1].spatial_smoothing_kernel != "none",
        )
        fuv = stack.enter_context(open_product(techniques[0][1].source_fuv_detector))
        frame_stop = min(stop_frame or fuv.shape[0], fuv.shape[0])
        frames = range(start_frame, frame_stop)
        for frame in tqdm(frames, desc=f"Orbit {orbit:04d} frames"):
            plot_frame(
                orbit,
                frame,
                fuv.time[frame],
                techniques,
                fuv,
                output,
                ratio_maximum,
            )


#%% Command line

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot detector counts and ratio for Product-2 smoothing variants."
    )
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument(
        "--product-directory",
        action="append",
        type=Path,
        help=(
            "Product-2 technique directory; repeat for each row. If omitted, "
            "IR_hardy and IR_hardy_smooth_* are discovered under --base."
        ),
    )
    parser.add_argument("--orbit", action="append", type=int)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--start-frame", type=int, default=0)
    parser.add_argument("--stop-frame", type=int)
    parser.add_argument("--ratio-maximum", type=float, default=RATIO_MAXIMUM)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    base = args.base.expanduser()
    output = args.output.expanduser()
    directories = discover_product_directories(base, args.product_directory)
    orbits = common_orbits(directories, args.orbit)
    if not orbits:
        raise FileNotFoundError("no common Product-2 orbit files found")

    print("Rows:")
    for directory in directories:
        print(f"  {directory}")
    print(f"Orbits: {orbits}")
    for orbit in orbits:
        plot_orbit(
            orbit,
            directories,
            output,
            args.start_frame,
            args.stop_frame,
            args.ratio_maximum,
        )
    print(f"Saved frame comparisons under {output}")


if __name__ == "__main__":
    main()
