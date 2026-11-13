"""Plot the compact FAST--IMAGE conjunction analysis dataset."""

#%% Imports and plot definitions

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from netCDF4 import Dataset
import numpy as np
from scipy.stats import rankdata


MLT_SLICES = (
    ("22-02", 22, 2), ("02-06", 2, 6), ("06-10", 6, 10),
    ("10-14", 10, 14), ("14-18", 14, 18), ("18-22", 18, 22),
)

ENERGY_RELATIONS = (
    ("full", "Full-spectrum mean energy", "pixel_earthward_mean_energy_full_kev"),
    ("500ev", ">=500 eV mean energy", "pixel_earthward_mean_energy_500ev_kev"),
    ("effective", "Robinson-effective energy", "pixel_robinson_effective_mean_energy_kev"),
)

FLUX_RELATIONS = (
    ("full", "Full-spectrum energy flux", "pixel_earthward_energy_flux_full_mw_m2"),
    ("500ev", ">=500 eV energy flux", "pixel_earthward_energy_flux_500ev_mw_m2"),
    ("effective", "Robinson-effective energy flux", "pixel_robinson_effective_energy_flux_mw_m2"),
)

CONDUCTANCE_RELATIONS = (
    ("hall", "Rees Hall conductance", "pixel_electron_hall_conductance_rees_s"),
    ("pedersen", "Rees Pedersen conductance", "pixel_electron_pedersen_conductance_rees_s"),
)

RATIO_EDGES = np.linspace(0, 150, 51)
WIC_EDGES = np.linspace(0, 2250, 61)
ENERGY_EDGES = np.linspace(0, 15, 61)
FLUX_EDGES = np.linspace(0, 25, 61)
HALL_EDGES = np.linspace(0, 60, 61)
PEDERSEN_EDGES = np.linspace(0, 30, 61)

SAMPLE_DISTRIBUTIONS = (
    ("sample_time_offset_seconds", np.linspace(-60, 60, 49), "Time offset [s]", "linear"),
    ("sample_detector_separation_deg", np.linspace(0, 1.3, 45), "Detector-centre separation [deg]", "linear"),
    ("sample_fast_altitude_km", np.linspace(400, 4100, 45), "FAST altitude [km]", "linear"),
    ("sample_footprint_mlt", np.linspace(0, 24, 49), "Footprint MLT [h]", "linear"),
    ("sample_footprint_mlat", np.linspace(-90, 90, 73), "Footprint MLAT [deg]", "linear"),
    ("sample_image_sza", np.linspace(0, 165, 56), "IMAGE SZA [deg]", "linear"),
    ("sample_earthward_energy_flux_all_valid_mw_m2", np.logspace(-4, 2, 61), "Full-spectrum energy flux [mW m$^{-2}$]", "log"),
    ("sample_earthward_energy_flux_centre_mw_m2", np.logspace(-4, 2, 61), ">=500 eV energy flux [mW m$^{-2}$]", "log"),
    ("sample_robinson_effective_energy_flux_mw_m2", np.logspace(-4, 2, 61), "Robinson-effective flux [mW m$^{-2}$]", "log"),
    ("sample_earthward_mean_energy_all_valid_kev", np.linspace(0, 15, 61), "Full-spectrum mean energy [keV]", "linear"),
    ("sample_earthward_mean_energy_centre_kev", np.linspace(0, 15, 61), ">=500 eV mean energy [keV]", "linear"),
    ("sample_robinson_effective_mean_energy_kev", np.linspace(0, 15, 61), "Robinson-effective energy [keV]", "linear"),
)

PIXEL_DISTRIBUTIONS = (
    ("pixel_time_offset_seconds", np.linspace(-60, 60, 49), "Time offset [s]", "linear"),
    ("pixel_detector_separation_deg", np.linspace(0, 1.3, 45), "Detector-centre separation [deg]", "linear"),
    ("pixel_fast_altitude_km", np.linspace(400, 4100, 45), "FAST altitude [km]", "linear"),
    ("pixel_footprint_mlt", np.linspace(0, 24, 49), "Footprint MLT [h]", "linear"),
    ("pixel_footprint_mlat", np.linspace(-90, 90, 73), "Footprint MLAT [deg]", "linear"),
    ("pixel_image_sza", np.linspace(0, 165, 56), "IMAGE SZA [deg]", "linear"),
    ("pixel_earthward_energy_flux_full_mw_m2", np.logspace(-4, 2, 61), "Full-spectrum energy flux [mW m$^{-2}$]", "log"),
    ("pixel_earthward_energy_flux_500ev_mw_m2", np.logspace(-4, 2, 61), ">=500 eV energy flux [mW m$^{-2}$]", "log"),
    ("pixel_robinson_effective_energy_flux_mw_m2", np.logspace(-4, 2, 61), "Robinson-effective flux [mW m$^{-2}$]", "log"),
    ("pixel_earthward_mean_energy_full_kev", np.linspace(0, 15, 61), "Full-spectrum mean energy [keV]", "linear"),
    ("pixel_earthward_mean_energy_500ev_kev", np.linspace(0, 15, 61), ">=500 eV mean energy [keV]", "linear"),
    ("pixel_robinson_effective_mean_energy_kev", np.linspace(0, 15, 61), "Robinson-effective energy [keV]", "linear"),
)


#%% Data and statistics

def values(variable, records=None):
    array = variable[:] if records is None else variable[records]
    if np.ma.isMaskedArray(array):
        if np.issubdtype(array.dtype, np.floating):
            array = array.filled(np.nan)
        else:
            array = array.filled(-1)
    return np.asarray(array)


def available(data, name):
    return name in data.variables


def correlation(x, y, rank=False):
    finite = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x, dtype=float)[finite]
    y = np.asarray(y, dtype=float)[finite]
    if x.size < 2:
        return np.nan
    if rank:
        x = rankdata(x)
        y = rankdata(y)
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def mlt_mask(mlt, start, stop):
    mlt = np.mod(np.asarray(mlt, dtype=float), 24)
    finite = np.isfinite(mlt)
    if start < stop:
        return finite & (mlt >= start) & (mlt < stop)
    return finite & ((mlt >= start) | (mlt < stop))


def pair_values(data, x_name, y_name, quality_name=None, sector=None,
                pixel_selected=None):
    x = values(data.variables[x_name]).astype(float)
    y = values(data.variables[y_name]).astype(float)
    selected = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if quality_name is not None and available(data, quality_name):
        selected &= values(data.variables[quality_name]).astype(bool)
    if pixel_selected is not None:
        selected &= pixel_selected
    if sector is not None:
        mlt = values(data.variables["pixel_footprint_mlt"])
        selected &= mlt_mask(mlt, sector[1], sector[2])
    return x[selected], y[selected]


def correlation_row(group, relation, sector, x_name, y_name, x, y):
    return {
        "group": group,
        "relation": relation,
        "mlt_sector": sector,
        "x_variable": x_name,
        "y_variable": y_name,
        "pair_count": int(x.size),
        "pearson_r": correlation(x, y),
        "spearman_r": correlation(x, y, rank=True),
        "log10_pearson_r": correlation(np.log10(x), np.log10(y)),
    }


def add_correlation(axis, x, y, logarithmic=False, fontsize=8):
    pearson = correlation(np.log10(x), np.log10(y)) if logarithmic else correlation(x, y)
    spearman = correlation(x, y, rank=True)
    label = "Pearson r (log)" if logarithmic else "Pearson r"
    axis.text(
        0.97, 0.97, f"{label} = {pearson:.2f}\nSpearman rho = {spearman:.2f}\nn = {x.size:,}",
        transform=axis.transAxes, ha="right", va="top", fontsize=fontsize,
        color="white", bbox={"facecolor": "black", "alpha": 0.55,
                             "edgecolor": "none", "pad": 2.5},
    )


def conditional_median(x, y, edges, minimum_count=20):
    index = np.digitize(x, edges) - 1
    median = np.full(edges.size - 1, np.nan)
    for bin_index in range(median.size):
        selected = index == bin_index
        if np.count_nonzero(selected) >= minimum_count:
            median[bin_index] = np.median(y[selected])
    return median


def write_csv(filename, rows):
    rows = list(rows)
    if not rows:
        return
    with Path(filename).open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def save_figure(figure, filename):
    figure.savefig(filename, dpi=180, bbox_inches="tight")
    plt.close(figure)


#%% Selection and distribution figures

def selection_rows(data, selected):
    year = values(data.variables["sample_year"])
    years = sorted(np.unique(year).tolist())
    rows = []
    cohorts = [(str(item), year == item) for item in years]
    cohorts.append(("all", np.ones(year.size, dtype=bool)))
    for cohort, in_cohort in cohorts:
        count = int(np.count_nonzero(in_cohort))
        pass_count = int(np.count_nonzero(selected & in_cohort))
        rows.append({
            "year": cohort, "criterion": "plot_filters",
            "label": "Plot filters", "sample_count": count,
            "pass_count": pass_count, "fail_count": count - pass_count,
            "pass_fraction": pass_count / count if count else np.nan,
        })
    return rows


def plot_selection(rows, filename):
    criteria = list(dict.fromkeys(row["criterion"] for row in rows))
    labels = [
        next(row["label"] for row in rows if row["criterion"] == criterion)
        for criterion in criteria
    ]
    years = sorted({row["year"] for row in rows if row["year"] != "all"})
    x = np.arange(len(criteria))
    width = 0.8 / max(1, len(years))
    figure, axis = plt.subplots(figsize=(13, 5.5))
    for index, year in enumerate(years):
        lookup = {
            row["criterion"]: row["pass_fraction"] for row in rows
            if row["year"] == year
        }
        axis.bar(
            x - 0.4 + width / 2 + index * width,
            [lookup.get(name, np.nan) for name in criteria],
            width=width, label=year,
        )
    axis.set_xticks(x, labels, rotation=25, ha="right")
    axis.set(ylabel="Fraction passing criterion", ylim=(0, 1.05),
             title="FAST--IMAGE selection criteria")
    axis.grid(axis="y", alpha=0.2)
    axis.legend(title="Year")
    save_figure(figure, filename)


def distribution_indices(count, maximum):
    if count <= maximum:
        return np.arange(count)
    return np.linspace(0, count - 1, maximum, dtype=np.int64)


def plot_distributions(data, settings, dimension, filename, title,
                       maximum_samples, separate_years=False,
                       selected=None):
    if selected is None:
        count = len(data.dimensions[dimension])
        records = distribution_indices(count, maximum_samples)
    else:
        available_records = np.flatnonzero(selected)
        count = available_records.size
        records = available_records[
            distribution_indices(count, maximum_samples)
        ]
    year_name = "sample_year" if dimension == "sample" else "pixel_year"
    year = values(data.variables[year_name], records)
    groups = sorted(np.unique(year)) if separate_years else [None]
    figure, axes = plt.subplots(4, 3, figsize=(13, 14))
    for axis, (name, bins, label, scale) in zip(axes.flat, settings):
        if not available(data, name):
            axis.text(0.5, 0.5, f"Unavailable: {name}", ha="center", va="center")
            axis.set_axis_off()
            continue
        array = values(data.variables[name], records).astype(float)
        for group in groups:
            selected = np.isfinite(array)
            if scale == "log":
                selected &= array > 0
            if group is not None:
                selected &= year == group
            if np.any(selected):
                axis.hist(array[selected], bins=bins, histtype="step",
                          density=True, linewidth=1.4,
                          label=str(group) if group is not None else None)
        axis.set(xlabel=label, ylabel="Density", xscale=scale)
        axis.grid(alpha=0.2)
    if separate_years:
        axes[0, 0].legend(title="Year")
    figure.suptitle(f"{title}; {count:,} {dimension}s")
    figure.tight_layout(rect=(0, 0, 1, 0.975))
    save_figure(figure, filename)


#%% Two-dimensional comparison figures

def plotted_relations(data, relations, x_name):
    return [relation for relation in relations
            if available(data, relation[2]) and available(data, x_name)]


def plot_fixed_x_relations(data, relations, x_name, x_edges, y_edges,
                           quality_name, filename, title, x_label, y_unit,
                           correlation_rows, group_name, by_mlt=False,
                           pixel_selected=None):
    relations = plotted_relations(data, relations, x_name)
    if not relations:
        return
    sectors = MLT_SLICES if by_mlt else (("all", 0, 24),)
    figure, axes = plt.subplots(
        len(relations), len(sectors),
        figsize=(3.5 * len(sectors) + 1.5, 3.1 * len(relations) + 1.2),
        squeeze=False, sharex=True, sharey=False, layout="constrained",
    )
    panel = {}
    maximum = 1
    for relation_name, label, y_name in relations:
        edges = y_edges[relation_name] if isinstance(y_edges, dict) else y_edges
        for sector in sectors:
            actual_sector = None if sector[0] == "all" else sector
            x, y = pair_values(
                data, x_name, y_name, quality_name, actual_sector,
                pixel_selected
            )
            count, _x, _y = np.histogram2d(x, y, bins=(x_edges, edges))
            panel[(relation_name, sector[0])] = (count, x, y, edges)
            maximum = max(maximum, int(np.max(count)))
            correlation_rows.append(correlation_row(
                group_name, relation_name, sector[0], x_name, y_name, x, y
            ))
    norm = LogNorm(vmin=1, vmax=max(2, maximum))
    x_centre = (x_edges[:-1] + x_edges[1:]) / 2
    plotted = None
    for row, (relation_name, label, _y_name) in enumerate(relations):
        for column, sector in enumerate(sectors):
            axis = axes[row, column]
            count, x, y, edges = panel[(relation_name, sector[0])]
            plotted = axis.pcolormesh(
                x_edges, edges, np.ma.masked_equal(count.T, 0),
                norm=norm, cmap="viridis", shading="flat",
            )
            axis.plot(x_centre, conditional_median(x, y, x_edges),
                      color="white", linewidth=1.4)
            add_correlation(axis, x, y, fontsize=7 if by_mlt else 8)
            axis.grid(alpha=0.12)
            if row == 0:
                axis.set_title(f"{sector[0]} MLT" if by_mlt else sector[0])
            if column == 0:
                axis.set_ylabel(f"{label}\n[{y_unit}]")
            if row == len(relations) - 1:
                axis.set_xlabel(x_label)
    figure.colorbar(plotted, ax=axes, label="IMAGE pixel count")
    figure.suptitle(title + "\nwhite line: conditional median (minimum 20 pixels)")
    save_figure(figure, filename)


def plot_image_fast(data, filename, correlation_rows, pixel_selected=None):
    definitions = []
    for name, label, fast_name in ENERGY_RELATIONS:
        definitions.append((f"energy_{name}", label, "pixel_image_E0", fast_name,
                            None, "keV"))
    for name, label, fast_name in FLUX_RELATIONS:
        definitions.append((f"flux_{name}", label, "pixel_image_Fe", fast_name,
                            None, "mW m-2"))
    definitions.extend((
        ("hall", "Hall conductance", "pixel_image_H",
         "pixel_electron_hall_conductance_rees_s", None, "S"),
        ("pedersen", "Pedersen conductance", "pixel_image_P",
         "pixel_electron_pedersen_conductance_rees_s", None, "S"),
    ))
    definitions = [item for item in definitions
                   if available(data, item[2]) and available(data, item[3])]
    if not definitions:
        return
    columns = 3
    rows = int(np.ceil(len(definitions) / columns))
    figure, axes = plt.subplots(rows, columns, figsize=(14, 4.3 * rows),
                                squeeze=False, layout="constrained")
    for axis, definition in zip(axes.flat, definitions):
        name, label, image_name, fast_name, quality_name, unit = definition
        image, fast = pair_values(
            data, image_name, fast_name, quality_name,
            pixel_selected=pixel_selected
        )
        if image.size:
            x = np.log10(image)
            y = np.log10(fast)
            axis.hexbin(x, y, gridsize=55, bins="log", mincnt=1, cmap="viridis")
            lower = min(np.percentile(x, 0.5), np.percentile(y, 0.5))
            upper = max(np.percentile(x, 99.5), np.percentile(y, 99.5))
            if lower == upper:
                lower -= 0.1
                upper += 0.1
            axis.plot([lower, upper], [lower, upper], color="white", linewidth=1)
            axis.set(xlim=(lower, upper), ylim=(lower, upper))
            add_correlation(axis, image, fast, logarithmic=True)
        axis.set(title=label, xlabel=f"log10 IMAGE [{unit}]",
                 ylabel=f"log10 FAST [{unit}]")
        axis.grid(alpha=0.12)
        correlation_rows.append(correlation_row(
            "image_fast", name, "all", image_name, fast_name, image, fast
        ))
    for axis in axes.flat[len(definitions):]:
        axis.set_axis_off()
    figure.suptitle("IMAGE products versus pixel-aggregated FAST")
    save_figure(figure, filename)


def plot_image_fast_by_mlt(data, relations, filename, title,
                           correlation_rows, group_name,
                           pixel_selected=None):
    relations = [item for item in relations
                 if available(data, item[2]) and available(data, item[3])]
    if not relations:
        return

    panel = {}
    maximum = 1
    for relation_name, _label, image_name, fast_name, edges, _unit in relations:
        for sector in MLT_SLICES:
            image, fast = pair_values(
                data, image_name, fast_name, sector=sector,
                pixel_selected=pixel_selected
            )
            count, _x, _y = np.histogram2d(
                image, fast, bins=(edges, edges)
            )
            panel[(relation_name, sector[0])] = count, image, fast
            maximum = max(maximum, int(np.max(count)))
            correlation_rows.append(correlation_row(
                group_name, relation_name, sector[0],
                image_name, fast_name, image, fast
            ))

    figure, axes = plt.subplots(
        len(relations), len(MLT_SLICES),
        figsize=(22.5, 3.2 * len(relations) + 1.2),
        squeeze=False, layout="constrained",
    )
    norm = LogNorm(vmin=1, vmax=max(2, maximum))
    plotted = None
    for row, (relation_name, label, _image_name, _fast_name,
              edges, unit) in enumerate(relations):
        centres = (edges[:-1] + edges[1:]) / 2
        for column, sector in enumerate(MLT_SLICES):
            axis = axes[row, column]
            count, image, fast = panel[(relation_name, sector[0])]
            plotted = axis.pcolormesh(
                edges, edges, np.ma.masked_equal(count.T, 0),
                norm=norm, cmap="viridis", shading="flat",
            )
            axis.plot(
                edges[[0, -1]], edges[[0, -1]],
                color="0.8", linestyle="--", linewidth=1,
            )
            axis.plot(
                centres, conditional_median(image, fast, edges),
                color="white", linewidth=1.4,
            )
            add_correlation(axis, image, fast, fontsize=6.5)
            axis.set(xlim=(edges[0], edges[-1]),
                     ylim=(edges[0], edges[-1]))
            axis.grid(alpha=0.12)
            if row == 0:
                axis.set_title(f"{sector[0]} MLT")
            if column == 0:
                axis.set_ylabel(f"FAST {label}\n[{unit}]")
            if row == len(relations) - 1:
                axis.set_xlabel(f"IMAGE [{unit}]")

    figure.colorbar(plotted, ax=axes, label="IMAGE pixel count")
    figure.suptitle(
        title + "\nwhite: conditional median; dashed: one-to-one"
    )
    save_figure(figure, filename)


def dombeck_diagnostic(data, filename, pixel_selected=None):
    name = "pixel_dombeck_intense_exposure_fraction"
    if not available(data, name):
        return
    fraction = values(data.variables[name]).astype(float)
    groups = (
        ("No intense exposure", np.isfinite(fraction) & (fraction == 0)),
        ("Any intense exposure", np.isfinite(fraction) & (fraction > 0)),
    )
    if pixel_selected is not None:
        groups = tuple(
            (label, selected & pixel_selected) for label, selected in groups
        )
    settings = (
        ("pixel_earthward_energy_flux_full_mw_m2", np.logspace(-4, 2, 61), "Full energy flux [mW m$^{-2}$]", "log"),
        ("pixel_earthward_mean_energy_full_kev", np.linspace(0, 15, 61), "Full mean energy [keV]", "linear"),
        ("pixel_earthward_mean_energy_500ev_kev", np.linspace(0, 15, 61), ">=500 eV mean energy [keV]", "linear"),
    )
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.4), layout="constrained")
    for axis, (variable, bins, label, scale) in zip(axes, settings):
        array = values(data.variables[variable]).astype(float)
        for group_label, selected in groups:
            valid = selected & np.isfinite(array)
            if scale == "log":
                valid &= array > 0
            axis.hist(array[valid], bins=bins, density=True, histtype="step",
                      linewidth=1.5, label=group_label)
        axis.set(xlabel=label, ylabel="Density", xscale=scale)
        axis.grid(alpha=0.2)
    axes[0].legend()
    figure.suptitle("Dombeck intensity diagnostic; not a quality selection")
    save_figure(figure, filename)


#%% Complete plotting workflow

def collection_rows(data, sample_selected, pixel_selected):
    sample_year = values(data.variables["sample_year"])
    pixel_year = values(data.variables["pixel_year"])
    rows = []
    for year in sorted(np.unique(sample_year)):
        sample = sample_year == year
        rows.append({
            "year": int(year),
            "sample_count": int(np.count_nonzero(sample)),
            "selected_sample_count": int(np.count_nonzero(
                sample & sample_selected
            )),
            "all_pixel_count": int(np.count_nonzero(pixel_year == year)),
            "pixel_count": int(np.count_nonzero(
                (pixel_year == year) & pixel_selected
            )),
        })
    return rows


def selection_masks(data, args):
    sample_selected = np.ones(len(data.dimensions["sample"]), dtype=bool)
    pixel_selected = np.ones(len(data.dimensions["pixel"]), dtype=bool)

    sample_valid = values(
        data.variables["sample_fast_measurement_valid"]
    ).astype(bool)
    sample_selected &= sample_valid
    pixel_valid_name = "pixel_fast_measurement_valid_exposure_fraction"
    if not available(data, pixel_valid_name):
        raise ValueError(
            "aggregated dataset lacks pixel FAST-valid exposure metadata; "
            "rebuild it with the current aggregate_conjunctions.py"
        )
    pixel_valid = values(data.variables[pixel_valid_name]).astype(float)
    pixel_selected &= np.isfinite(pixel_valid) & (pixel_valid > 0)

    annotation = args.annotation_status
    if annotation != "all":
        code = {"accepted": 1, "rejected": 0, "unannotated": -1}[annotation]
        sample_status = values(
            data.variables["sample_annotation_status"]
        ).astype(np.int8)
        pixel_status = values(
            data.variables["pixel_annotation_status"]
        ).astype(np.int8)
        sample_selected &= sample_status == code
        pixel_selected &= pixel_status == code

    population = args.fast_population
    if population != "all-valid":
        sample_available = values(
            data.variables["sample_dombeck_available"]
        ).astype(bool)
        sample_intense = values(
            data.variables["sample_dombeck_intense_precipitation"]
        ).astype(bool)
        sample_selected &= sample_available

        available_fraction = values(data.variables[
            "pixel_dombeck_available_exposure_fraction"
        ]).astype(float)
        intense_fraction = values(data.variables[
            "pixel_dombeck_intense_exposure_fraction"
        ]).astype(float)
        pixel_selected &= np.isfinite(available_fraction)
        pixel_selected &= available_fraction > 0
        pixel_selected &= np.isfinite(intense_fraction)
        if population == "dombeck-intense":
            sample_selected &= sample_intense
            pixel_selected &= intense_fraction > 0
        else:
            sample_selected &= ~sample_intense
            pixel_selected &= intense_fraction == 0

    filters = (
        ("maximum_detector_separation_deg",
         "sample_detector_separation_deg",
         "pixel_maximum_detector_separation_deg", np.less_equal, False),
        ("maximum_altitude_km",
         "sample_fast_altitude_km",
         "pixel_maximum_fast_altitude_km", np.less_equal, False),
        ("minimum_sza_deg", "sample_image_sza",
         "pixel_image_sza", np.greater_equal, False),
        ("maximum_dza_deg", "sample_image_dza",
         "pixel_image_dza", np.less_equal, False),
        ("maximum_absolute_time_offset_seconds",
         "sample_time_offset_seconds",
         "pixel_maximum_absolute_time_offset_seconds", np.less_equal, True),
    )
    for argument, sample_name, pixel_name, comparison, absolute in filters:
        limit = getattr(args, argument)
        if limit is None:
            continue
        if not available(data, sample_name) or not available(data, pixel_name):
            raise ValueError(
                f"aggregated dataset lacks {sample_name} or {pixel_name}"
            )
        sample_array = values(data.variables[sample_name]).astype(float)
        pixel_array = values(data.variables[pixel_name]).astype(float)
        if absolute:
            sample_array = np.abs(sample_array)
            pixel_array = np.abs(pixel_array)
        sample_selected &= np.isfinite(sample_array)
        sample_selected &= comparison(sample_array, limit)
        pixel_selected &= np.isfinite(pixel_array)
        pixel_selected &= comparison(pixel_array, limit)
    return sample_selected, pixel_selected


def build_plots(args):
    output = Path(args.output)
    if output.exists() and not output.is_dir():
        raise NotADirectoryError(output)
    if output.exists() and any(output.iterdir()) and not args.overwrite:
        raise FileExistsError(f"{output} is not empty; pass --overwrite")
    output.mkdir(parents=True, exist_ok=True)

    with Dataset(args.input) as data:
        if "processing_complete" not in data.ncattrs():
            raise ValueError("compact conjunction dataset has no completion flag")
        if data.getncattr("processing_complete") != "true":
            raise ValueError("compact conjunction dataset is incomplete")
        if "sample" not in data.dimensions or "pixel" not in data.dimensions:
            raise ValueError("input is not a compact conjunction analysis dataset")

        sample_selected, pixel_selected = selection_masks(data, args)
        selections = selection_rows(data, sample_selected)
        write_csv(output / "selection_summary.csv", selections)
        write_csv(
            output / "collection_summary.csv",
            collection_rows(data, sample_selected, pixel_selected),
        )
        plot_selection(selections, output / "selection_overview.png")
        plot_distributions(
            data, SAMPLE_DISTRIBUTIONS, "sample",
            output / "distribution_before.png",
            "Selected FAST--IMAGE sample distribution",
            args.maximum_distribution_samples, selected=sample_selected,
        )
        plot_distributions(
            data, PIXEL_DISTRIBUTIONS, "pixel",
            output / "distribution_after.png",
            "Selected pixel-aggregated FAST--IMAGE distribution",
            args.maximum_distribution_samples, selected=pixel_selected,
        )
        plot_distributions(
            data, PIXEL_DISTRIBUTIONS, "pixel",
            output / "coverage_overview.png",
            "Selected pixel distribution by year",
            args.maximum_distribution_samples, separate_years=True,
            selected=pixel_selected,
        )

        correlations = []
        plot_image_fast(
            data, output / "image_fast_comparison.png", correlations,
            pixel_selected,
        )
        energy_pairs = tuple(
            (name, label, "pixel_image_E0", fast_name, ENERGY_EDGES, "keV")
            for name, label, fast_name in ENERGY_RELATIONS
        )
        plot_image_fast_by_mlt(
            data, energy_pairs, output / "energy_comparison_by_mlt.png",
            "IMAGE-derived energy versus pixel-aggregated FAST energy by MLT",
            correlations, "image_fast_energy_mlt", pixel_selected,
        )
        conductance_pairs = (
            (
                "hall", "Hall conductance", "pixel_image_H",
                "pixel_electron_hall_conductance_rees_s", HALL_EDGES, "S",
            ),
            (
                "pedersen", "Pedersen conductance", "pixel_image_P",
                "pixel_electron_pedersen_conductance_rees_s",
                PEDERSEN_EDGES, "S",
            ),
        )
        plot_image_fast_by_mlt(
            data, conductance_pairs,
            output / "conductance_comparison_by_mlt.png",
            "IMAGE-derived versus FAST Rees conductance by MLT",
            correlations, "image_fast_conductance_mlt", pixel_selected,
        )
        plot_fixed_x_relations(
            data, ENERGY_RELATIONS, "pixel_image_R", RATIO_EDGES, ENERGY_EDGES,
            None, output / "ratio_energy_comparison.png",
            "IMAGE WIC/SI13 ratio versus FAST mean energy",
            "IMAGE corrected WIC/SI13 ratio", "keV", correlations,
            "ratio_energy", pixel_selected=pixel_selected,
        )
        plot_fixed_x_relations(
            data, ENERGY_RELATIONS, "pixel_image_R", RATIO_EDGES, ENERGY_EDGES,
            None, output / "ratio_energy_by_mlt.png",
            "IMAGE WIC/SI13 ratio versus FAST mean energy by MLT",
            "IMAGE corrected WIC/SI13 ratio", "keV", correlations,
            "ratio_energy", by_mlt=True, pixel_selected=pixel_selected,
        )

        groups = (
            ("energy", ENERGY_RELATIONS, ENERGY_EDGES, "keV"),
            ("flux", FLUX_RELATIONS, FLUX_EDGES, "mW m-2"),
            ("conductance", CONDUCTANCE_RELATIONS,
             {"hall": HALL_EDGES, "pedersen": PEDERSEN_EDGES}, "S"),
        )
        for group_name, relations, edges, unit in groups:
            plot_fixed_x_relations(
                data, relations, "pixel_image_wic_corrected", WIC_EDGES,
                edges, None, output / f"wic_{group_name}_comparison.png",
                f"Corrected IMAGE WIC versus FAST {group_name}",
                "IMAGE corrected WIC [counts]", unit, correlations,
                f"wic_{group_name}", pixel_selected=pixel_selected,
            )
            plot_fixed_x_relations(
                data, relations, "pixel_image_wic_corrected", WIC_EDGES,
                edges, None, output / f"wic_{group_name}_by_mlt.png",
                f"Corrected IMAGE WIC versus FAST {group_name} by MLT",
                "IMAGE corrected WIC [counts]", unit, correlations,
                f"wic_{group_name}", by_mlt=True,
                pixel_selected=pixel_selected,
            )
        dombeck_diagnostic(
            data, output / "dombeck_diagnostic.png", pixel_selected
        )
        write_csv(output / "correlation_summary.csv", correlations)

    print(f"Saved FAST--IMAGE analysis figures and tables: {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--maximum-distribution-samples", type=int,
                        default=250_000)
    parser.add_argument(
        "--annotation-status",
        choices=("all", "accepted", "rejected", "unannotated"),
        default="all",
    )
    parser.add_argument(
        "--fast-population",
        choices=("all-valid", "dombeck-intense", "dombeck-weak"),
        default="all-valid",
    )
    parser.add_argument("--maximum-detector-separation-deg", type=float)
    parser.add_argument("--maximum-altitude-km", type=float)
    parser.add_argument("--minimum-sza-deg", type=float)
    parser.add_argument("--maximum-dza-deg", type=float)
    parser.add_argument("--maximum-absolute-time-offset-seconds", type=float)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(args.input)
    if args.maximum_distribution_samples < 1:
        raise ValueError("--maximum-distribution-samples must be positive")
    build_plots(args)


if __name__ == "__main__":
    main()
