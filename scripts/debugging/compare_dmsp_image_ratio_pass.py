"""Compare DMSP SSJ and IMAGE-ratio energies over one complete polar pass."""

#%% Imports

import argparse
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from icreader import load as icload
import numpy as np
import xarray as xr

from icphysics.image import wic_to_s13


#%% Read the processed products

def read_image(filename):
    """Read the image-ratio precipitation product and its stored grid."""

    data = icload(filename)
    if data.method != "image_ratio":
        raise ValueError("The IMAGE file must use method='image_ratio'")

    return {
        "time": np.asarray(data.time, dtype="datetime64[ns]"),
        "ratio": np.asarray(data.R, dtype=float),
        "energy": np.asarray(data.E0, dtype=float),
        "xi": np.asarray(data.grid.xi, dtype=float),
        "eta": np.asarray(data.grid.eta, dtype=float),
        "projection": data.grid.projection,
    }


def read_dmsp(filename):
    """Read the complete mapped DMSP pass from the reduced product."""

    with xr.open_dataset(filename) as data:
        return {
            "time": np.asarray(data.time.values, dtype="datetime64[ns]"),
            "mlat": np.asarray(data.footprint_qd_lat.values, dtype=float),
            "mlt": np.asarray(data.footprint_mlt.values, dtype=float),
            "energy": np.asarray(data.ELE_AVG_ENERGY.values, dtype=float) / 1000,
        }


#%% Sample IMAGE along the DMSP pass

def nearest_grid_cells(xi, eta, track_xi, track_eta):
    """Find the nearest regular Cubed-Sphere cell for each track point."""

    xi_index = np.argmin(np.abs(xi[0, :, None] - track_xi[None, :]), axis=0)
    eta_index = np.argmin(np.abs(eta[:, 0, None] - track_eta[None, :]), axis=0)
    return eta_index, xi_index


def sample_image(image, dmsp, maximum_time_gap_seconds):
    """Sample the nearest IMAGE frame and grid cell along the DMSP track."""

    frame_index = np.argmin(
        np.abs(dmsp["time"][:, None] - image["time"][None, :]),
        axis=1,
    )
    time_gap = np.abs(dmsp["time"] - image["time"][frame_index])

    track_xi, track_eta = image["projection"].geo2cube(
        dmsp["mlt"] * 15,
        dmsp["mlat"],
    )
    eta_index, xi_index = nearest_grid_cells(
        image["xi"], image["eta"], track_xi, track_eta
    )

    inside_grid = (
        (track_xi >= image["xi"].min())
        & (track_xi <= image["xi"].max())
        & (track_eta >= image["eta"].min())
        & (track_eta <= image["eta"].max())
    )
    close_in_time = time_gap <= np.timedelta64(maximum_time_gap_seconds, "s")
    valid = inside_grid & close_in_time

    ratio = image["ratio"][frame_index, eta_index, xi_index]
    energy = image["energy"][frame_index, eta_index, xi_index]
    ratio[~valid] = np.nan
    energy[~valid] = np.nan

    return ratio, energy, valid


def moving_average(values, width):
    """Calculate a centred, NaN-aware moving average."""

    finite = np.isfinite(values)
    total = np.convolve(np.where(finite, values, 0), np.ones(width), mode="same")
    count = np.convolve(finite.astype(float), np.ones(width), mode="same")

    smooth = np.full(values.shape, np.nan)
    use = count >= width // 2 + 1
    smooth[use] = total[use] / count[use]
    return smooth


#%% Figure

def save_figure(dmsp, ratio, image_energy, output, smoothing_seconds):
    """Plot the complete DMSP pass and the processed IMAGE comparison."""

    dmsp_smooth = moving_average(dmsp["energy"], smoothing_seconds)
    figure, axes = plt.subplots(
        3, 1, figsize=(11, 8), sharex=True,
        gridspec_kw={"height_ratios": [2, 1, 1]},
        constrained_layout=True,
    )

    axes[0].plot(dmsp["time"], dmsp["energy"], color="0.75", linewidth=0.7, label="DMSP SSJ, 1 s")
    axes[0].plot(dmsp["time"], dmsp_smooth, color="black", linewidth=1.8, label=f"DMSP, {smoothing_seconds}-s mean")
    plotted_energy = np.minimum(image_energy, 15)
    axes[0].plot(dmsp["time"], plotted_energy, color="crimson", linewidth=1.5, label="IMAGE ratio method, capped at 15 keV")
    axes[0].set_ylabel("Mean energy [keV]")
    axes[0].set_ylim(0, 15)
    axes[0].legend(frameon=False, loc="upper left")

    axes[1].plot(dmsp["time"], ratio, color="tab:purple", linewidth=1.3)
    axes[1].axhline(np.max(wic_to_s13), color="black", linestyle="--", linewidth=1, label="Frey table maximum")
    axes[1].set_ylabel("WIC/SI13")
    axes[1].set_yscale("log")
    axes[1].legend(frameon=False, loc="upper left")

    axes[2].plot(dmsp["time"], dmsp["mlat"], color="tab:blue", linewidth=1.5)
    axes[2].axhline(40, color="0.5", linestyle=":", linewidth=1)
    axes[2].set_ylabel("QD latitude [deg]")
    axes[2].set_xlabel("UTC")

    for axis in axes:
        axis.grid(alpha=0.2)
    time_locator = mdates.AutoDateLocator(minticks=6, maxticks=10)
    axes[2].xaxis.set_major_locator(time_locator)
    axes[2].xaxis.set_major_formatter(mdates.ConciseDateFormatter(time_locator))

    figure.suptitle("DMSP F15 and IMAGE WIC/SI13 over the complete polar pass")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output.with_suffix(".png"), dpi=200)
    figure.savefig(output.with_suffix(".pdf"))
    plt.close(figure)


#%% Command line

def parse_args():
    repository = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--image-file", type=Path,
        default=Path("/home/bing/dtu_server/IMAGE_FUV/=precipitation_IR_P2/or_0968.nc"),
    )
    parser.add_argument(
        "--dmsp-file", type=Path,
        default=repository / "example_data/dmsp/dmsp-f15_ssj_20011021-20011022_north.nc",
    )
    parser.add_argument(
        "--output", type=Path,
        default=repository / "figures/debugging/coumans_2004/dmsp_image_ratio_complete_pass",
    )
    parser.add_argument("--maximum-time-gap-seconds", type=int, default=75)
    parser.add_argument("--smoothing-seconds", type=int, default=9)
    return parser.parse_args()


def main():
    args = parse_args()
    image = read_image(args.image_file)
    dmsp = read_dmsp(args.dmsp_file)
    ratio, energy, valid = sample_image(
        image, dmsp, args.maximum_time_gap_seconds
    )
    save_figure(dmsp, ratio, energy, args.output, args.smoothing_seconds)

    print(f"DMSP interval: {dmsp['time'][0]} to {dmsp['time'][-1]}")
    print(f"DMSP samples: {dmsp['time'].size}")
    print(f"IMAGE comparison samples: {valid.sum()}")
    print(f"Saved {args.output.with_suffix('.png')}")
    print(f"Saved {args.output.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
