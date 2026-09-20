"""Calculate precipitation on the WIC detector geometry."""

#%% Imports

from pathlib import Path

import numpy as np
from icphysics import (
    PROTON_RESPONSE_ENERGY_RANGE,
    hardy_ion_precipitation,
    precipitation_from_ratio,
    proton_correct_images,
)
from netCDF4 import Dataset, date2num, num2date

from .kp import load_gfz_kp, match_gfz_kp
from .fuvdetector import (
    PREPROCESSING_LABEL,
    SCHEMA_VERSION as FUV_SCHEMA_VERSION,
    SOURCE_TIME_DECODING,
    source_identity,
)


#%% Product configuration

SCHEMA_VERSION = 2
PRECIPITATION_METHOD = "image_ratio"
PROTON_ENERGY_MODELS = ("hardy", "constant")
PROTON_FLUX_SOURCE = "SI12"
PROTON_OPERATION_ORDER = (
    "SI12 counts are coregistered onto WIC pixels in fuv_detector, then "
    "converted to proton flux on the WIC detector geometry"
)
COUNT_UNCERTAINTY_METHOD = (
    "square roots of Product-1 detector-count variances are supplied to "
    "icPhysics. SI variances already include independent-pixel propagation "
    "through coregistration; covariance and background-model uncertainty are "
    "not included"
)
HARDY_COORDINATE_NOTE = (
    "Hardy corrected geomagnetic coordinates approximated by Product-1 "
    "Modified Apex latitude and MLT at 130 km"
)


FRAME_FIELDS = (
    "wic_counts", "si12_counts", "si13_counts",
    "wic_variance", "si12_variance", "si13_variance",
    "wic_quality_weight", "si12_quality_weight", "si13_quality_weight",
    "wic_coverage", "si12_coverage", "si13_coverage",
    "wic_valid", "si12_valid", "si13_valid",
    "si12_source_count", "si13_source_count",
    "glat", "glon", "mlat", "mlon", "mlt", "sza", "dza",
)
TIME_FIELDS = (
    "time", "wic_source_time", "si12_source_time", "si13_source_time",
    "wic_source_index", "si12_source_index", "si13_source_index", "ssalon",
    "wic_frame_quality", "si12_frame_quality", "si13_frame_quality",
)


#%% Product-1 loading

def _as_array(variable, dtype=float):
    values = variable[:]
    if np.ma.isMaskedArray(values):
        if np.issubdtype(np.dtype(dtype), np.integer):
            fill_value = -1
        elif np.issubdtype(np.dtype(dtype), np.bool_):
            fill_value = False
        else:
            fill_value = np.nan
        values = values.filled(fill_value)
    return np.asarray(values, dtype=dtype)


def _read_time(nc, name):
    variable = nc.variables[name]
    encoded = variable[:]
    missing = np.ma.getmaskarray(encoded)
    decoded = np.full(encoded.shape, None, dtype=object)
    if np.any(~missing):
        decoded[~missing] = num2date(
            np.asarray(encoded)[~missing],
            variable.units,
            variable.calendar,
            only_use_cftime_datetimes=False,
        )
    return decoded


def load_fuv_detector(filename):
    """Load the Product-1 fields needed by detector precipitation."""

    filename = Path(filename)
    with Dataset(filename) as nc:
        if (
            nc.product_type != "fuv_detector"
            or nc.representation != "detector"
            or int(nc.schema_version) != FUV_SCHEMA_VERSION
            or nc.preprocessing_label != PREPROCESSING_LABEL
            or getattr(nc, "source_time_decoding", None) != SOURCE_TIME_DECODING
        ):
            raise ValueError(
                f"{filename} is not a supported schema-{FUV_SCHEMA_VERSION} "
                f"{PREPROCESSING_LABEL} fuv_detector product"
            )

        product = {
            "source_file": str(filename),
            "schema_version": int(nc.schema_version),
            "preprocessing_label": nc.preprocessing_label,
            "source_time_decoding": nc.source_time_decoding,
            "source_software_version": nc.software_version,
            "coordinate_system": nc.coordinate_system,
            "reference_height_km": float(nc.reference_height_km),
            "shape": (
                len(nc.dimensions["time"]),
                len(nc.dimensions["row"]),
                len(nc.dimensions["column"]),
            ),
            "detector_row": _as_array(nc.variables["detector_row"], int),
            "detector_column": _as_array(
                nc.variables["detector_column"], int
            ),
        }
        for name in TIME_FIELDS:
            if name.endswith("_time") or name == "time":
                product[name] = _read_time(nc, name)
            elif name.endswith("_index") or name.endswith("_frame_quality"):
                product[name] = _as_array(nc.variables[name], int)
            else:
                product[name] = _as_array(nc.variables[name])
        for name in FRAME_FIELDS:
            dtype = bool if name.endswith("_valid") else float
            product[name] = _as_array(nc.variables[name], dtype)

    product.update(source_identity(filename))
    return product


#%% Proton energy

def make_detector_proton_energy(
    kp,
    mlt,
    mlat,
    model,
    constant_energy,
    constant_uncertainty,
):
    """Evaluate proton energy on each time-dependent WIC detector pixel."""

    if model not in PROTON_ENERGY_MODELS:
        raise ValueError(f"proton_energy_model must be one of {PROTON_ENERGY_MODELS}")
    if constant_uncertainty < 0:
        raise ValueError("proton_energy_uncertainty must not be negative")

    shape = np.asarray(mlt).shape
    if model == "constant":
        if not np.isfinite(constant_energy) or constant_energy <= 0:
            raise ValueError("proton_energy must be a positive finite value")
        Ep_model = np.full(shape, constant_energy, dtype=float)
        dEp = np.full(shape, constant_uncertainty, dtype=float)
        uncertainty_method = "constant supplied value"
        coordinate_note = "not applicable to constant proton energy"
    else:
        Ep_model = np.full(shape, np.nan)
        for frame, kp_value in enumerate(np.asarray(kp, dtype=float)):
            valid = (
                np.isfinite(kp_value)
                & np.isfinite(mlt[frame])
                & np.isfinite(mlat[frame])
            )
            if not np.any(valid):
                continue
            hardy = hardy_ion_precipitation(
                kp_value,
                mlt[frame][valid],
                mlat[frame][valid],
            )
            Ep_model[frame][valid] = hardy["mean_energy"]
        dEp = np.where(np.isfinite(Ep_model), 0.0, np.nan)
        uncertainty_method = "not modelled by Hardy et al. (1991)"
        coordinate_note = HARDY_COORDINATE_NOTE

    lower, upper = PROTON_RESPONSE_ENERGY_RANGE
    Ep = np.clip(Ep_model, lower, upper)
    clipping_flag = np.isfinite(Ep_model) & (Ep != Ep_model)

    return Ep_model, Ep, dEp, clipping_flag, uncertainty_method, coordinate_note


#%% Detector precipitation

class PrecipitationDetector:
    """Apply the image-ratio precipitation method to one fuv_detector orbit."""

    def __init__(
        self,
        fuv_detector,
        *,
        kp_series=None,
        proton_energy_model="hardy",
        proton_energy=2.0,
        proton_energy_uncertainty=0.0,
        software_version="unrecorded experimental worktree",
    ):
        fuv = load_fuv_detector(fuv_detector)
        if kp_series is None:
            kp_series = load_gfz_kp()

        self.product_type = "precipitation_detector"
        self.representation = "detector"
        self.schema_version = SCHEMA_VERSION
        self.method = PRECIPITATION_METHOD
        self.proton_flux_source = PROTON_FLUX_SOURCE
        self.proton_energy_model = proton_energy_model
        self.proton_energy_constant = float(proton_energy)
        self.proton_energy_uncertainty_constant = float(
            proton_energy_uncertainty
        )
        self.proton_operation_order = PROTON_OPERATION_ORDER
        self.count_uncertainty_method = COUNT_UNCERTAINTY_METHOD
        self.software_version = str(software_version)
        self.source_fuv_detector = fuv["source_file"]
        self.source_fuv_detector_sha256 = fuv["source_sha256"]
        self.source_fuv_detector_size_bytes = fuv["source_size_bytes"]
        self.source_fuv_detector_mtime_ns = fuv["source_mtime_ns"]
        self.source_fuv_detector_schema_version = fuv["schema_version"]
        self.source_preprocessing_label = fuv["preprocessing_label"]
        self.source_fuv_detector_time_decoding = fuv["source_time_decoding"]
        self.source_software_version = fuv["source_software_version"]
        self.coordinate_system = fuv["coordinate_system"]
        self.reference_height_km = fuv["reference_height_km"]
        self.shape = fuv["shape"]

        for name in TIME_FIELDS + FRAME_FIELDS + (
            "detector_row", "detector_column",
        ):
            setattr(self, name, np.asarray(fuv[name]))

        # 1. Match authoritative Kp to every retained WIC frame.
        matched_kp = match_gfz_kp(self.time, kp_series)
        self.kp = matched_kp["kp"]
        self.kp_interval_start = matched_kp["interval_start"]
        self.kp_provenance = dict(kp_series["provenance"])

        # 2. Evaluate proton energy on the time-dependent detector geometry.
        proton_energy = make_detector_proton_energy(
            self.kp,
            self.mlt,
            self.mlat,
            self.proton_energy_model,
            self.proton_energy_constant,
            self.proton_energy_uncertainty_constant,
        )
        (
            self.Ep_model,
            self.Ep,
            self.dEp,
            self.Ep_clipping_flag,
            self.proton_energy_uncertainty_method,
            self.proton_energy_coordinate_note,
        ) = proton_energy

        # 3. Require all three image-ratio channels at each detector pixel.
        input_valid = (
            self.wic_valid
            & self.si12_valid
            & self.si13_valid
            & np.isfinite(self.Ep)
        )
        wic = np.where(input_valid, self.wic_counts, np.nan)
        si12 = np.where(input_valid, self.si12_counts, np.nan)
        si13 = np.where(input_valid, self.si13_counts, np.nan)

        # 4. Infer proton flux from mapped SI12, then correct WIC and SI13.
        dwic = np.where(input_valid, np.sqrt(self.wic_variance), np.nan)
        dsi12 = np.where(input_valid, np.sqrt(self.si12_variance), np.nan)
        dsi13 = np.where(input_valid, np.sqrt(self.si13_variance), np.nan)
        with np.errstate(divide="ignore", invalid="ignore"):
            corrected = proton_correct_images(
                wic=wic,
                dwic=dwic,
                si12=si12,
                dsi12=dsi12,
                si13=si13,
                dsi13=dsi13,
                proton_energy=self.Ep,
                proton_energy_uncertainty=self.dEp,
            )

        for name, values in corrected.items():
            setattr(self, name, np.asarray(values, dtype=float))

        # 5. Calculate image-ratio electron energy and energy flux.
        # Keep central values when a propagated uncertainty becomes undefined.
        dwic_for_ratio = np.where(
            np.isfinite(self.dwic_corrected), self.dwic_corrected, 0.0
        )
        dsi13_for_ratio = np.where(
            np.isfinite(self.dsi13_corrected), self.dsi13_corrected, 0.0
        )
        ratio = precipitation_from_ratio(
            self.wic_corrected,
            dwic_for_ratio,
            self.si13_corrected,
            dsi13_for_ratio,
        )
        for name, values in ratio.items():
            setattr(self, name, np.asarray(values, dtype=float))

        # 6. Retain input support separately from successful method output.
        self.method_valid = input_valid & np.isfinite(self.E0) & np.isfinite(self.Fe)
        self.method_quality_weight = np.where(
            self.method_valid,
            self.wic_quality_weight
            * self.si12_quality_weight
            * self.si13_quality_weight,
            np.nan,
        )
        uncertain = ~(
            np.isfinite(self.dwic_corrected)
            & np.isfinite(self.dsi13_corrected)
            & self.method_valid
        )
        for name in ("dR", "dE0", "dFe"):
            values = getattr(self, name).copy()
            values[uncertain] = np.nan
            setattr(self, name, values)

        # Raw counts remain in the referenced Product 1 and are not duplicated
        # in Product 2 after the precipitation calculation is complete.
        del self.wic_counts, self.si12_counts, self.si13_counts

    #%% NetCDF output

    @staticmethod
    def _write_time(nc, name, values, units):
        values = np.asarray(values)
        if np.issubdtype(values.dtype, np.datetime64):
            values = values.astype("datetime64[ms]").astype(object)
        encoded = np.full(values.shape, np.nan)
        valid = np.asarray([value is not None for value in values])
        if np.any(valid):
            encoded[valid] = date2num(
                values[valid].tolist(), units, calendar="standard"
            )
        variable = nc.createVariable(name, "f8", ("time",), fill_value=np.nan)
        variable[:] = encoded
        variable.units = units
        variable.calendar = "standard"
        variable.time_zone = "UTC"

    def to_nc(self, filename):
        """Write one self-contained detector-space precipitation orbit."""

        with Dataset(filename, "w", format="NETCDF4") as nc:
            nc.createDimension("time", self.shape[0])
            nc.createDimension("row", self.shape[1])
            nc.createDimension("column", self.shape[2])

            nc.product_type = self.product_type
            nc.representation = self.representation
            nc.schema_version = self.schema_version
            nc.method = self.method
            nc.proton_flux_source = self.proton_flux_source
            nc.proton_energy_model = self.proton_energy_model
            nc.proton_operation_order = self.proton_operation_order
            nc.proton_energy_uncertainty_method = (
                self.proton_energy_uncertainty_method
            )
            nc.proton_energy_coordinate_note = self.proton_energy_coordinate_note
            nc.count_uncertainty_method = self.count_uncertainty_method
            nc.proton_response_energy_min = PROTON_RESPONSE_ENERGY_RANGE[0]
            nc.proton_response_energy_max = PROTON_RESPONSE_ENERGY_RANGE[1]
            nc.source_fuv_detector = self.source_fuv_detector
            nc.source_fuv_detector_sha256 = self.source_fuv_detector_sha256
            nc.source_fuv_detector_size_bytes = np.int64(
                self.source_fuv_detector_size_bytes
            )
            nc.source_fuv_detector_mtime_ns = np.int64(
                self.source_fuv_detector_mtime_ns
            )
            nc.source_fuv_detector_schema_version = (
                self.source_fuv_detector_schema_version
            )
            nc.source_preprocessing_label = self.source_preprocessing_label
            nc.source_fuv_detector_time_decoding = (
                self.source_fuv_detector_time_decoding
            )
            nc.source_software_version = self.source_software_version
            nc.coordinate_system = self.coordinate_system
            nc.reference_height_km = self.reference_height_km
            nc.software_version = self.software_version
            if self.proton_energy_model == "constant":
                nc.proton_energy_constant = self.proton_energy_constant
                nc.proton_energy_uncertainty_constant = (
                    self.proton_energy_uncertainty_constant
                )
            for name, value in self.kp_provenance.items():
                nc.setncattr(f"kp_{name}", value)

            time_units = "seconds since 2000-01-01 00:00:00"
            time_fields = (
                ("time", "time"),
                ("wic_source_time", "wic_source_time"),
                ("si12_source_time", "si12_source_time"),
                ("si13_source_time", "si13_source_time"),
                ("Kp_interval_start", "kp_interval_start"),
            )
            for output_name, attribute_name in time_fields:
                self._write_time(
                    nc,
                    output_name,
                    getattr(self, attribute_name),
                    time_units,
                )

            for name in (
                "wic_source_index", "si12_source_index", "si13_source_index"
            ):
                variable = nc.createVariable(name, "i4", ("time",))
                variable[:] = getattr(self, name)

            for name in (
                "wic_frame_quality", "si12_frame_quality", "si13_frame_quality"
            ):
                variable = nc.createVariable(
                    name, "i1", ("time",), fill_value=-1
                )
                variable[:] = getattr(self, name)
                variable.flag_values = np.asarray([0, 1, 2], dtype=np.int8)
                variable.flag_meanings = "rejected usable science_ready"

            for name in ("detector_row", "detector_column"):
                dimension = "row" if name == "detector_row" else "column"
                variable = nc.createVariable(name, "i4", (dimension,))
                variable[:] = getattr(self, name)

            time_fields = {
                "Kp": (self.kp, "1"),
                "ssalon": (self.ssalon, "degrees"),
            }
            for name, (values, units) in time_fields.items():
                variable = nc.createVariable(name, "f4", ("time",))
                variable[:] = values
                variable.units = units

            dimensions = ("time", "row", "column")
            fields = {
                "glat": (self.glat, "degrees_north"),
                "glon": (self.glon, "degrees_east"),
                "mlat": (self.mlat, "degrees"),
                "mlon": (self.mlon, "degrees"),
                "mlt": (self.mlt, "hours"),
                "sza": (self.sza, "degrees"),
                "dza": (self.dza, "degrees"),
                "wic_quality_weight": (self.wic_quality_weight, "1"),
                "si12_quality_weight": (self.si12_quality_weight, "1"),
                "si13_quality_weight": (self.si13_quality_weight, "1"),
                "method_quality_weight": (self.method_quality_weight, "1"),
                "wic_coverage": (self.wic_coverage, "1"),
                "si12_coverage": (self.si12_coverage, "1"),
                "si13_coverage": (self.si13_coverage, "1"),
                "Ep_model": (self.Ep_model, "keV"),
                "Ep": (self.Ep, "keV"),
                "dEp": (self.dEp, "keV"),
                "Fp": (self.Fp, "mW m-2"),
                "dFp": (self.dFp, "mW m-2"),
                "wic_corrected": (self.wic_corrected, "counts"),
                "dwic_corrected": (self.dwic_corrected, "counts"),
                "si13_corrected": (self.si13_corrected, "counts"),
                "dsi13_corrected": (self.dsi13_corrected, "counts"),
                "R": (self.R, "1"),
                "dR": (self.dR, "1"),
                "E0": (self.E0, "keV"),
                "dE0": (self.dE0, "keV"),
                "Fe": (self.Fe, "mW m-2"),
                "dFe": (self.dFe, "mW m-2"),
                "varE0Fe": (self.varE0Fe, "keV mW m-2"),
            }
            for name, (values, units) in fields.items():
                dtype = "f8" if name in ("glat", "glon", "mlat", "mlon", "mlt") else "f4"
                variable = nc.createVariable(name, dtype, dimensions, zlib=True)
                variable[:] = values
                variable.units = units

            for name in ("wic_valid", "si12_valid", "si13_valid", "method_valid"):
                variable = nc.createVariable(name, "i1", dimensions, zlib=True)
                variable[:] = getattr(self, name).astype(np.int8)

            clipped = nc.createVariable(
                "Ep_clipping_flag", "i1", dimensions, zlib=True
            )
            clipped[:] = self.Ep_clipping_flag.astype(np.int8)

            for name in ("si12_source_count", "si13_source_count"):
                variable = nc.createVariable(name, "i4", dimensions, zlib=True)
                variable[:] = getattr(self, name)
