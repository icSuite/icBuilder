"""Construct method-specific precipitation products from native binned images."""

#%% Imports

from copy import deepcopy
from datetime import timedelta
from pathlib import Path

import numpy as np
from icreader import load as icload
from netCDF4 import Dataset, date2num
from icphysics import (
    PROTON_RESPONSE_ENERGY_RANGE,
    hardy_ion_precipitation,
    precipitation_from_ratio,
    precipitation_from_hardy,
    precipitation_from_zhang_paxton,
    proton_correct_images,
)

from .grids import grid_mlt
from .precipitationdetector import ELECTRON_METADATA_ATTRIBUTES, hardy_electron_metadata
from .kp import load_gfz_kp, match_gfz_kp


#%% Helpers

def load_binned_image(image):
    """Load a binned NetCDF filename, or return an existing image object."""

    if isinstance(image, (str, Path)):
        return icload(image)
    return image

#%% Time matching

TIME_TOLERANCE = timedelta(seconds=2)
TIME_MATCH_RULE = (
    "Image ratio matches WIC, SI12, and SI13 one-to-one. Zhang-Paxton matches "
    "WIC and SI12 one-to-one, as does Hardy, then attaches SI13 without reducing frame "
    "support. Every accepted sensor group spans at most 2 seconds; the latest "
    "required-sensor timestamp is the product time. Missing SI13 index is -1."
)

def match_sensor_times(*sensor_times):
    """Match two or three sensor time series using the legacy orbit rule."""

    sorted_times = [sorted((time, index) for index, time in enumerate(times)) for times in sensor_times]
    pointers = [0] * len(sorted_times)
    common_times = []
    source_indices = []

    while all(pointer < len(times) for pointer, times in zip(pointers, sorted_times)):
        current = [times[pointer] for pointer, times in zip(pointers, sorted_times)]
        current_times = [item[0] for item in current]

        earliest = min(current_times)
        latest = max(current_times)

        if latest - earliest <= TIME_TOLERANCE:
            common_times.append(latest)
            source_indices.append([item[1] for item in current])
            pointers = [pointer + 1 for pointer in pointers]
        else:
            # None of the tied earliest frames can match the current later
            # frame, so discard all of them together.
            for sensor, time in enumerate(current_times):
                if time == earliest:
                    pointers[sensor] += 1

    return (np.asarray(common_times, dtype=object), 
            np.asarray(source_indices, dtype=int).reshape(-1, len(sensor_times)))


#%% Spatial interpolation

REGRID_METHOD = (
    "Direct copy when source and target grids are identical; otherwise "
    "bilinear interpolation between regular Cubed-Sphere cell centres and "
    "nearest source cell at the physical grid boundary"
)
REGRID_UNCERTAINTY = (
    "All four surrounding source cells must be finite. Target variance is "
    "sum(a_i^2 sigma_i^2) using bilinear weights a_i and assuming independent "
    "source-cell errors. Boundary targets inherit the nearest source-cell "
    "uncertainty. Covariance between target cells is not stored."
)

def make_regrid_mapping(source_grid, target_grid):
    """Calculate the fixed bilinear mapping between two regular grids."""

    source_x_grid = np.asarray(source_grid.xi)
    source_y_grid = np.asarray(source_grid.eta)
    source_x = source_x_grid[0, :]
    source_y = source_y_grid[:, 0]

    if not np.allclose(source_x_grid, source_x[None, :]):
        raise ValueError("source-grid xi must be constant down each column")
    if not np.allclose(source_y_grid, source_y[:, None]):
        raise ValueError("source-grid eta must be constant across each row")
    if np.any(np.diff(source_x) <= 0) or np.any(np.diff(source_y) <= 0):
        raise ValueError("source-grid xi and eta must increase")

    target_x = np.asarray(target_grid.xi).ravel()
    target_y = np.asarray(target_grid.eta).ravel()

    # Find the four source cells surrounding each target cell.
    x_upper = np.clip(np.searchsorted(source_x, target_x), 1, len(source_x) - 1)
    y_upper = np.clip(np.searchsorted(source_y, target_y), 1, len(source_y) - 1)
    x_lower = x_upper - 1
    y_lower = y_upper - 1

    inside_centres = (
        (target_x >= source_x[0]) & (target_x <= source_x[-1])
        & (target_y >= source_y[0]) & (target_y <= source_y[-1])
    )
    target_indices = np.flatnonzero(inside_centres)

    x_fraction = ((target_x[inside_centres] - source_x[x_lower[inside_centres]]) /
                  (source_x[x_upper[inside_centres]] - source_x[x_lower[inside_centres]]))
    
    y_fraction = ((target_y[inside_centres] - source_y[y_lower[inside_centres]]) /
                  (source_y[y_upper[inside_centres]] - source_y[y_lower[inside_centres]]))

    weights = np.column_stack([(1 - x_fraction) * (1 - y_fraction), 
                               x_fraction * (1 - y_fraction), 
                               (1 - x_fraction) * y_fraction, 
                               x_fraction * y_fraction])
    
    rows = np.column_stack([y_lower[inside_centres],
                            y_lower[inside_centres],
                            y_upper[inside_centres],
                            y_upper[inside_centres]])
    
    columns = np.column_stack([x_lower[inside_centres],
                               x_upper[inside_centres],
                               x_lower[inside_centres],
                               x_upper[inside_centres]])

    # Cell-centre interpolation does not reach the outer half of the edge
    # cells. Those locations are still inside the physical source grid and
    # inherit the nearest edge-cell value. This must not fill internal gaps.
    x_edges = (source_x[0] - np.diff(source_x)[0] / 2,
               source_x[-1] + np.diff(source_x)[-1] / 2)
    y_edges = (source_y[0] - np.diff(source_y)[0] / 2,
               source_y[-1] + np.diff(source_y)[-1] / 2)
    inside_edges = (
        (target_x >= x_edges[0]) & (target_x <= x_edges[1])
        & (target_y >= y_edges[0]) & (target_y <= y_edges[1])
    )
    boundary = inside_edges & ~inside_centres
    boundary_indices = np.flatnonzero(boundary)

    nearest_columns = np.abs(
        target_x[boundary, None] - source_x[None, :]
    ).argmin(axis=1)
    nearest_rows = np.abs(
        target_y[boundary, None] - source_y[None, :]
    ).argmin(axis=1)

    return {"source_shape": source_grid.shape,
            "target_shape": target_grid.shape,
            "target_indices": target_indices,
            "rows": rows,
            "columns": columns,
            "weights": weights,
            "boundary_indices": boundary_indices,
            "nearest_rows": nearest_rows,
            "nearest_columns": nearest_columns}


def regrid_to_target(values, mapping, propagate_uncertainty=False):
    """Apply a bilinear grid mapping to one source array."""

    values = np.asarray(values, dtype=float)
    if values.ndim != 3 or values.shape[1:] != mapping["source_shape"]:
        raise ValueError("source values do not match the source grid")

    target_shape = (values.shape[0], *mapping["target_shape"])
    target_values = np.full(target_shape, np.nan)

    target_indices = mapping["target_indices"]
    rows = mapping["rows"]
    columns = mapping["columns"]
    weights = mapping["weights"]

    for frame in range(values.shape[0]):
        surrounding_values = values[frame, rows, columns]
        valid = np.all(np.isfinite(surrounding_values), axis=1)

        selected_weights = weights[valid]
        selected_targets = target_indices[valid]

        if propagate_uncertainty:
            interpolated_values = np.sqrt(np.sum(selected_weights**2 * surrounding_values[valid]**2,axis=1))
        else:
            interpolated_values = np.sum(selected_weights * surrounding_values[valid], axis=1)

        target_values[frame].flat[selected_targets] = interpolated_values

        # Fill only the geometric boundary. A missing nearest source cell
        # remains NaN, and internal missing-data gaps are never filled here.
        boundary_values = values[
            frame,
            mapping["nearest_rows"],
            mapping["nearest_columns"],
        ]
        target_values[frame].flat[mapping["boundary_indices"]] = boundary_values

    return target_values


def grids_are_identical(source_grid, target_grid):
    """Return True when two grids describe the same cell centres."""

    if source_grid.shape != target_grid.shape:
        return False

    return (
        np.allclose(source_grid.xi, target_grid.xi)
        and np.allclose(source_grid.eta, target_grid.eta)
        and np.allclose(source_grid.lat, target_grid.lat)
        and np.allclose(source_grid.lon, target_grid.lon)
    )


def prepare_si_sensor(sensor, indices, target_grid):
    """Place one SI sensor on the Product-2 time axis and WIC grid."""

    shape = (len(indices), *target_grid.shape)
    values = np.full(shape, np.nan)
    sigma = np.full(shape, np.nan)
    weight = np.full(shape, np.nan)

    matched = indices >= 0
    if not np.any(matched):
        return values, sigma, weight

    sensor_indices = indices[matched]

    # Do not let four-neighbour interpolation enlarge holes when the SI data
    # are already on exactly the requested grid. Preserve every value, NaN,
    # uncertainty, and weight cell-for-cell.
    if grids_are_identical(sensor.grid, target_grid):
        values[matched] = np.asarray(sensor.mu)[sensor_indices]
        sigma[matched] = np.asarray(sensor.sigma)[sensor_indices]
        weight[matched] = np.asarray(sensor.w)[sensor_indices]
        return values, sigma, weight

    mapping = make_regrid_mapping(sensor.grid, target_grid)
    values[matched] = regrid_to_target(sensor.mu[sensor_indices], mapping)
    sigma[matched] = regrid_to_target(sensor.sigma[sensor_indices], mapping, propagate_uncertainty=True)
    weight[matched] = regrid_to_target(sensor.w[sensor_indices], mapping)

    return values, sigma, weight


def match_optional_si13(wic_times, si12_times, si13_times):
    """Attach SI13 without allowing it to remove matched WIC/SI12 frames."""

    si13_times = sorted(
        (time, index) for index, time in enumerate(si13_times)
    )
    si13_indices = np.full(len(wic_times), -1, dtype=int)
    pointer = 0

    for frame, (wic_time, si12_time) in enumerate(zip(wic_times, si12_times)):
        earliest = min(wic_time, si12_time)
        latest = max(wic_time, si12_time)
        lower = latest - TIME_TOLERANCE
        upper = earliest + TIME_TOLERANCE

        while pointer < len(si13_times) and si13_times[pointer][0] < lower:
            pointer += 1

        if pointer < len(si13_times) and si13_times[pointer][0] <= upper:
            si13_indices[frame] = si13_times[pointer][1]
            pointer += 1

    return si13_indices


def prepare_image_ratio(wic, si12, si13, target_grid):
    """Prepare the three simultaneously matched image-ratio inputs."""

    time, indices = match_sensor_times(wic.time, si12.time, si13.time)
    source_indices = {"wic": indices[:, 0], "si12": indices[:, 1], "si13": indices[:, 2]}
    return prepare_observations(time, source_indices, wic, si12, si13, target_grid)


def prepare_zhang_paxton(wic, si12, si13, target_grid):
    """Prepare WIC/SI12 support and attach SI13 where it is available."""

    time, indices = match_sensor_times(wic.time, si12.time)
    source_indices = {"wic": indices[:, 0], "si12": indices[:, 1], "si13": np.full(len(time), -1, dtype=int)}

    if si13 is not None:
        wic_times = np.asarray(wic.time)[source_indices["wic"]]
        si12_times = np.asarray(si12.time)[source_indices["si12"]]
        source_indices["si13"] = match_optional_si13(wic_times, si12_times, si13.time)

    return prepare_observations(time, source_indices, wic, si12, si13, target_grid)


def prepare_observations(time, source_indices, wic, si12, si13, target_grid):
    """Collect matched WIC and regridded SI arrays in one dictionary."""

    wic_index = source_indices["wic"]
    prepared = {"time": time, 
                "source_indices": source_indices,
                "wic": np.asarray(wic.mu)[wic_index].copy(),
                "dwic": np.asarray(wic.sigma)[wic_index].copy(),
                "wic_weight": np.asarray(wic.w)[wic_index].copy()}

    si12_index = source_indices["si12"]
    prepared["si12"], prepared["dsi12"], prepared["si12_weight"] = (prepare_si_sensor(si12, si12_index, target_grid))

    if si13 is None:
        shape = (len(time), *target_grid.shape)
        prepared["si13"] = np.full(shape, np.nan)
        prepared["dsi13"] = np.full(shape, np.nan)
        prepared["si13_weight"] = np.full(shape, np.nan)
    else:
        si13_index = source_indices["si13"]
        prepared["si13"], prepared["dsi13"], prepared["si13_weight"] = (prepare_si_sensor(si13, si13_index, target_grid))

    return prepared


#%% Proton energy

PROTON_ENERGY_MODELS = ("hardy", "constant")
PROTON_FLUX_SOURCE = "SI12"
HARDY_COORDINATE_NOTE = (
    "Hardy corrected geomagnetic coordinates approximated by the Product-2 "
    "Modified Apex grid at 130 km"
)


def make_proton_energy(kp, grid, shape, model, constant_energy, constant_uncertainty):
    """Return raw and response-clipped proton mean energy on the image grid."""

    if model not in PROTON_ENERGY_MODELS:
        raise ValueError(f"proton_energy_model must be one of {PROTON_ENERGY_MODELS}")
    if constant_uncertainty < 0:
        raise ValueError("proton_energy_uncertainty must not be negative")

    if model == "constant":
        if not np.isfinite(constant_energy) or constant_energy <= 0:
            raise ValueError("proton_energy must be a positive finite value")
        Ep_model = np.full(shape, constant_energy, dtype=float)
        dEp = np.full(shape, constant_uncertainty, dtype=float)
        uncertainty_method = "constant supplied value"
        coordinate_note = "not applicable to constant proton energy"
    else:
        Ep_model = np.full(shape, np.nan)
        mlt = grid_mlt(grid)
        mlat = np.asarray(grid.lat, dtype=float)

        # Kp is constant over each three-hour GFZ interval. Evaluate each
        # distinct Hardy map once, then copy it into the matching frames.
        for kp_value in np.unique(kp[np.isfinite(kp)]):
            frames = kp == kp_value
            hardy = hardy_ion_precipitation(kp_value, mlt, mlat)
            Ep_model[frames] = hardy["mean_energy"]

        # Hardy supplies no predictive uncertainty for mean proton energy.
        # Zero here means "not modelled", not known without uncertainty.
        dEp = np.where(np.isfinite(Ep_model), 0.0, np.nan)
        uncertainty_method = "not modelled by Hardy et al. (1991)"
        coordinate_note = HARDY_COORDINATE_NOTE

    lower, upper = PROTON_RESPONSE_ENERGY_RANGE
    Ep = np.clip(Ep_model, lower, upper)
    Ep_clipping_flag = np.isfinite(Ep_model) & (Ep != Ep_model)

    return {
        "Ep_model": Ep_model,
        "Ep": Ep,
        "dEp": dEp,
        "Ep_clipping_flag": Ep_clipping_flag,
        "uncertainty_method": uncertainty_method,
        "coordinate_note": coordinate_note,
    }


#%% Precipitation product

class PrecipitationImage:
    """Prepare binned sensor data and call an injected precipitation model."""

    def __init__(
        self,
        wic,
        si12,
        method,
        physics_function=None,
        *,
        si13=None,
        kp_series=None,
        proton_energy_model="hardy",
        proton_energy=2.0,
        proton_energy_uncertainty=0.0,
        proton_flux_source=PROTON_FLUX_SOURCE,
        source_products=None,
        physics_provenance=None):
        
        if method not in ("image_ratio", "zhang_paxton", "hardy"):
            raise ValueError("method must be image_ratio, zhang_paxton, or hardy")
        if proton_flux_source != PROTON_FLUX_SOURCE:
            raise ValueError("proton_flux_source must currently be 'SI12'")
        if physics_function is None:
            physics_function = {
                "image_ratio": precipitation_from_ratio,
                "hardy": precipitation_from_hardy,
                "zhang_paxton": precipitation_from_zhang_paxton,
            }[method]

        input_files = {}
        for name, image in (("wic", wic), ("si12", si12), ("si13", si13)):
            if isinstance(image, (str, Path)):
                input_files[name] = str(image)

        wic = load_binned_image(wic)
        si12 = load_binned_image(si12)
        if si13 is not None:
            si13 = load_binned_image(si13)

        if wic.sensor != "WIC" or si12.sensor != "SI12":
            raise ValueError("wic and si12 inputs must contain WIC and SI12 data")
        if method == "image_ratio" and (si13 is None or si13.sensor != "SI13"):
            raise ValueError("image_ratio requires a SI13 BinnedImage")
        if si13 is not None and si13.sensor != "SI13":
            raise ValueError("si13 input must contain SI13 data")
        if kp_series is None:
            kp_series = load_gfz_kp()

        self.method = method
        if method == "hardy":
            for name, value in hardy_electron_metadata(HARDY_COORDINATE_NOTE).items():
                setattr(self, name, value)
        self.grid = deepcopy(wic.grid)
        self.proton_flux_source = str(proton_flux_source)
        self.proton_energy_model = str(proton_energy_model)
        self.proton_response_energy_min, self.proton_response_energy_max = PROTON_RESPONSE_ENERGY_RANGE
        self.proton_energy_constant = float(proton_energy)
        self.proton_energy_uncertainty_constant = float(proton_energy_uncertainty)
        self.source_products = input_files
        self.source_products.update(source_products or {})

        provenance = {"module": physics_function.__module__, 
                      "function": getattr(physics_function, "__name__", type(physics_function).__name__)}  # What is this?
        provenance.update(physics_provenance or {})
        self.physics_provenance = provenance

        # Each method states explicitly which sensors define its time support.
        if method == "image_ratio":
            prepared = prepare_image_ratio(wic, si12, si13, self.grid)
        else:
            prepared = prepare_zhang_paxton(wic, si12, si13, self.grid)

        self.time = prepared.pop("time")
        if self.time.size == 0:
            raise ValueError("the required sensors have no matching frames")
        self.source_indices = prepared.pop("source_indices")
        self.shape = (self.time.size, *self.grid.shape)
        self.ssalon = np.asarray(wic.ssalon)[self.source_indices["wic"]].copy()
        for name, values in prepared.items():
            setattr(self, name, values)

        # Weight only the sensors used by this precipitation method. SI13 is
        # diagnostic in the Zhang--Paxton product and must not reduce support.
        if method == "image_ratio":
            self.w = self.wic_weight * self.si12_weight * self.si13_weight
        else:
            self.w = self.wic_weight * self.si12_weight

        # Kp is assigned after the final method-specific time match.
        matched_kp = match_gfz_kp(self.time, kp_series)
        self.kp = matched_kp["kp"]
        self.kp_interval_start = matched_kp["interval_start"]
        self.kp_provenance = dict(kp_series["provenance"])

        proton_energy_fields = make_proton_energy(
            self.kp,
            self.grid,
            self.shape,
            self.proton_energy_model,
            self.proton_energy_constant,
            self.proton_energy_uncertainty_constant,
        )
        self.proton_energy_uncertainty_method = proton_energy_fields.pop(
            "uncertainty_method"
        )
        self.proton_energy_coordinate_note = proton_energy_fields.pop(
            "coordinate_note"
        )
        for name, values in proton_energy_fields.items():
            setattr(self, name, np.asarray(values).copy())

        corrected = proton_correct_images(wic=self.wic, dwic=self.dwic,
                                          si12=self.si12, dsi12=self.dsi12,
                                          si13=self.si13, dsi13=self.dsi13,
                                          proton_energy=self.Ep,
                                          proton_energy_uncertainty=self.dEp,
                                          **({"uncertainty_mode": "measurement"} if method == "hardy" else {}))
        
        for name, values in corrected.items():
            setattr(self, name, np.asarray(values, dtype=float).copy())

        inputs = {"wic_corrected": self.wic_corrected,"dwic_corrected": self.dwic_corrected}
        
        if method == "image_ratio":
            inputs.update(si13_corrected=self.si13_corrected, dsi13_corrected=self.dsi13_corrected)
        elif method == "hardy":
            inputs.update(kp=self.kp[:, None, None], mlt=grid_mlt(self.grid), mlat=np.asarray(self.grid.lat))
        else:
            inputs.update(kp=self.kp, mlt=grid_mlt(self.grid))

        result = physics_function(**inputs)
        if not isinstance(result, dict):
            raise TypeError("physics_function must return a dictionary")
        for forbidden in ("counts", "P", "H", "dP", "dH", "wic_corrected", "si13_corrected"):
            if forbidden in result:
                raise ValueError(f"physics_function must not return {forbidden}")

        required = ["E0", "dE0", "Fe", "dFe", "varE0Fe"]

        for name in required:
            if name not in result:
                raise ValueError(f"physics_function did not return {name}")
            values = np.asarray(result[name], dtype=float)
            if values.shape != self.shape:
                raise ValueError(f"{name} must have shape {self.shape}")
            setattr(self, name, values.copy())

        for name in ("R", "dR"):
            if name in result:
                values = np.asarray(result[name], dtype=float)
                if values.shape != self.shape:
                    raise ValueError(f"{name} must have shape {self.shape}")
                setattr(self, name, values.copy())

        if method == "hardy":
            self.E0_response_range_flag = np.asarray(result["E0_response_range_flag"], dtype=np.int8).copy()
            if self.E0_response_range_flag.shape != self.shape:
                raise ValueError(f"E0_response_range_flag must have shape {self.shape}")

        self.sensor_provenance = {}
        for name, sensor in (("wic", wic), ("si12", si12), ("si13", si13)):
            if sensor is None:
                continue
            self.sensor_provenance[name] = {"image_correction": sensor.correction or "raw", 
                                            "los_correction": bool(sensor.los_correction)}

    #%% NetCDF output

    def to_nc(self, filename):
        """Save the prepared observations and precipitation estimates."""

        with Dataset(filename, "w", format="NETCDF4") as nc:
            nc.createDimension("time", self.shape[0])
            nc.createDimension("dim1", self.shape[1])
            nc.createDimension("dim2", self.shape[2])

            nc.product_type = "precipitation"
            nc.schema_version = 2
            nc.method = self.method
            if self.method == "hardy":
                for name in ELECTRON_METADATA_ATTRIBUTES:
                    nc.setncattr(name, getattr(self, name))
            nc.proton_flux_source = self.proton_flux_source
            nc.proton_energy_model = self.proton_energy_model
            nc.proton_energy_uncertainty_method = self.proton_energy_uncertainty_method
            nc.proton_energy_coordinate_note = self.proton_energy_coordinate_note
            nc.proton_response_energy_min = PROTON_RESPONSE_ENERGY_RANGE[0]
            nc.proton_response_energy_max = PROTON_RESPONSE_ENERGY_RANGE[1]
            if self.proton_energy_model == "constant":
                nc.proton_energy_constant = self.proton_energy_constant
                nc.proton_energy_uncertainty_constant = self.proton_energy_uncertainty_constant
            nc.time_match_tolerance_seconds = TIME_TOLERANCE.total_seconds()
            nc.time_match_rule = TIME_MATCH_RULE
            nc.regrid_method = REGRID_METHOD
            nc.regrid_uncertainty = REGRID_UNCERTAINTY

            for name, value in self.physics_provenance.items():
                nc.setncattr(f"physics_{name}", value)
            for name, value in self.kp_provenance.items():
                nc.setncattr(f"kp_{name}", value)
            for sensor, provenance in self.sensor_provenance.items():
                nc.setncattr(f"{sensor}_image_correction", provenance["image_correction"])
                nc.setncattr(f"{sensor}_los_correction", np.int8(provenance["los_correction"]))
            for sensor, path in self.source_products.items():
                nc.setncattr(f"source_{sensor}", str(path))

            time_units = "seconds since 2000-01-01 00:00:00"
            time = nc.createVariable("time", "f8", ("time",))
            time[:] = date2num(self.time.tolist(), time_units, calendar="standard")
            time.units = time_units
            time.calendar = "standard"
            time.time_zone = "UTC"

            for sensor, indices in self.source_indices.items():
                variable = nc.createVariable(f"{sensor}_source_index", "i4", ("time",))
                variable[:] = indices

            kp = nc.createVariable("Kp", "f4", ("time",))
            kp[:] = self.kp
            kp.units = "1"

            kp_start = nc.createVariable("Kp_interval_start", "f8", ("time",))
            kp_start[:] = date2num(
                self.kp_interval_start.astype("datetime64[ms]").astype(object).tolist(),
                time_units,
                calendar="standard",
            )
            kp_start.units = time_units
            kp_start.calendar = "standard"

            ssalon = nc.createVariable("ssalon", "f4", ("time",))
            ssalon[:] = self.ssalon
            ssalon.units = "degrees"

            fields = {
                "wic": (self.wic, "counts"),
                "dwic": (self.dwic, "counts"),
                "si12": (self.si12, "counts"),
                "dsi12": (self.dsi12, "counts"),
                "si13": (self.si13, "counts"),
                "dsi13": (self.dsi13, "counts"),
                "wic_weight": (self.wic_weight, "1"),
                "si12_weight": (self.si12_weight, "1"),
                "si13_weight": (self.si13_weight, "1"),
                "w": (self.w, "1"),
                "wic_corrected": (self.wic_corrected, "counts"),
                "dwic_corrected": (self.dwic_corrected, "counts"),
                "si13_corrected": (self.si13_corrected, "counts"),
                "dsi13_corrected": (self.dsi13_corrected, "counts"),
                "Ep_model": (self.Ep_model, "keV"),
                "Ep": (self.Ep, "keV"),
                "dEp": (self.dEp, "keV"),
                "Fp": (self.Fp, "mW m-2"),
                "dFp": (self.dFp, "mW m-2"),
                "E0": (self.E0, "keV"),
                "dE0": (self.dE0, "keV"),
                "Fe": (self.Fe, "mW m-2"),
                "dFe": (self.dFe, "mW m-2"),
                "varE0Fe": (self.varE0Fe, "keV mW m-2"),
            }
            for name in ("R", "dR"):
                if hasattr(self, name):
                    fields[name] = (getattr(self, name), "1")

            for name, (data, units) in fields.items():
                variable = nc.createVariable(
                    name, "f4", ("time", "dim1", "dim2"), zlib=True
                )
                variable[:] = data
                variable.units = units

            if self.method == "hardy":
                flag = nc.createVariable("E0_response_range_flag", "i1", ("time", "dim1", "dim2"), zlib=True)
                flag[:] = self.E0_response_range_flag
                flag.flag_values = np.asarray([-1, 0, 1, 2], dtype=np.int8)
                flag.flag_meanings = "unavailable in_range below above"

            clipped = nc.createVariable(
                "Ep_clipping_flag", "i1", ("time", "dim1", "dim2"), zlib=True
            )
            clipped[:] = self.Ep_clipping_flag.astype(np.int8)
            clipped.long_name = "proton energy clipped to camera-response table range"

            grid = nc.createGroup("grid")
            grid.position = self.grid.projection.position.astype(float)
            grid.orientation = self.grid.projection.orientation
            grid.L = self.grid.L
            grid.W = self.grid.W
            grid.Lres = self.grid.Lres
            grid.Wres = self.grid.Wres
            grid.R = self.grid.R

            coordinates = {
                "xi": (self.grid.xi, "radians"),
                "eta": (self.grid.eta, "radians"),
                "mlat": (self.grid.lat, "degrees"),
                "mlt": (grid_mlt(self.grid), "hours"),
            }
            for name, (data, units) in coordinates.items():
                variable = grid.createVariable(name, "f8", ("dim1", "dim2"), zlib=True)
                variable[:] = data
                variable.units = units
