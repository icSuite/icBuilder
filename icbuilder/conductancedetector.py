"""Calculate conductance on the WIC detector geometry."""

#%% Imports

from pathlib import Path

import numpy as np
from icphysics import robinson_conductance
from netCDF4 import Dataset, date2num, num2date

from .fuvdetector import SOURCE_TIME_DECODING, source_identity
from .precipitationdetector import (
    PRECIPITATION_METHOD,
    SCHEMA_VERSION as PRECIPITATION_SCHEMA_VERSION,
)


#%% Product configuration

SCHEMA_VERSION = 1
CONDUCTANCE_MODEL = "robinson"
CONDUCTANCE_UNCERTAINTY_METHOD = (
    "first-order Robinson propagation of Product-2 dE0, dFe, and varE0Fe; "
    "at Fe=0, dP and dH are one-sided excursions from zero to dFe"
)

TIME_FIELDS = (
    "time", "wic_source_time", "si12_source_time", "si13_source_time",
    "Kp_interval_start",
)
TIME_VALUE_FIELDS = (
    "wic_source_index", "si12_source_index", "si13_source_index",
    "wic_frame_quality", "si12_frame_quality", "si13_frame_quality",
    "Kp", "ssalon",
)
GEOMETRY_FIELDS = (
    "glat", "glon", "mlat", "mlon", "mlt", "sza", "dza",
)
PRECIPITATION_FIELDS = (
    "Ep_model", "Ep", "dEp", "Fp", "dFp",
    "E0", "dE0", "Fe", "dFe", "varE0Fe",
)


#%% Product-2 loading

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


def load_precipitation_detector(filename):
    """Load the Product-2 state needed by detector conductance."""

    filename = Path(filename)
    with Dataset(filename) as nc:
        if (
            nc.product_type != "precipitation_detector"
            or nc.representation != "detector"
            or int(nc.schema_version) != PRECIPITATION_SCHEMA_VERSION
            or nc.method != PRECIPITATION_METHOD
            or getattr(nc, "source_fuv_detector_time_decoding", None)
            != SOURCE_TIME_DECODING
        ):
            raise ValueError(
                f"{filename} is not a supported schema-"
                f"{PRECIPITATION_SCHEMA_VERSION} {PRECIPITATION_METHOD} "
                "precipitation_detector product"
            )

        shape = (
            len(nc.dimensions["time"]),
            len(nc.dimensions["row"]),
            len(nc.dimensions["column"]),
        )
        if any(length == 0 for length in shape):
            raise ValueError("precipitation_detector dimensions must be non-empty")
        if nc.variables["detector_row"].shape != (shape[1],):
            raise ValueError(f"detector_row must have shape {(shape[1],)}")
        if nc.variables["detector_column"].shape != (shape[2],):
            raise ValueError(
                f"detector_column must have shape {(shape[2],)}"
            )
        product = {
            "source_file": str(filename),
            "schema_version": int(nc.schema_version),
            "method": nc.method,
            "proton_flux_source": nc.proton_flux_source,
            "proton_energy_model": nc.proton_energy_model,
            "proton_energy_uncertainty_method": (
                nc.proton_energy_uncertainty_method
            ),
            "proton_energy_coordinate_note": nc.proton_energy_coordinate_note,
            "proton_response_energy_min": float(
                nc.proton_response_energy_min
            ),
            "proton_response_energy_max": float(
                nc.proton_response_energy_max
            ),
            "proton_operation_order": nc.proton_operation_order,
            "count_uncertainty_method": nc.count_uncertainty_method,
            "coordinate_system": nc.coordinate_system,
            "reference_height_km": float(nc.reference_height_km),
            "source_software_version": nc.software_version,
            "source_fuv_detector": nc.source_fuv_detector,
            "source_fuv_detector_sha256": nc.source_fuv_detector_sha256,
            "source_preprocessing_label": nc.source_preprocessing_label,
            "source_fuv_detector_time_decoding": (
                nc.source_fuv_detector_time_decoding
            ),
            "shape": shape,
            "detector_row": _as_array(nc.variables["detector_row"], int),
            "detector_column": _as_array(
                nc.variables["detector_column"], int
            ),
            "kp_provenance": {
                name[3:]: nc.getncattr(name)
                for name in nc.ncattrs()
                if name.startswith("kp_")
            },
        }
        if product["proton_energy_model"] == "constant":
            product["proton_energy_constant"] = float(
                nc.proton_energy_constant
            )
            product["proton_energy_uncertainty_constant"] = float(
                nc.proton_energy_uncertainty_constant
            )

        for name in TIME_FIELDS:
            if nc.variables[name].shape != (shape[0],):
                raise ValueError(f"{name} must have shape {(shape[0],)}")
            product[name] = _read_time(nc, name)

        for name in TIME_VALUE_FIELDS:
            if nc.variables[name].shape != (shape[0],):
                raise ValueError(f"{name} must have shape {(shape[0],)}")
            dtype = int if name.endswith("_index") or name.endswith(
                "_frame_quality"
            ) else float
            product[name] = _as_array(nc.variables[name], dtype)

        for name in GEOMETRY_FIELDS + PRECIPITATION_FIELDS + (
            "method_quality_weight",
        ):
            if nc.variables[name].shape != shape:
                raise ValueError(f"{name} must have shape {shape}")
            product[name] = _as_array(nc.variables[name])

        for name in ("method_valid", "Ep_clipping_flag"):
            if nc.variables[name].shape != shape:
                raise ValueError(f"{name} must have shape {shape}")
            product[name] = _as_array(nc.variables[name], bool)

    product.update(source_identity(filename))
    return product


#%% Detector conductance

class ConductanceDetector:
    """Apply Robinson conductance to one detector Product-2 orbit."""

    def __init__(
        self,
        precipitation_detector,
        *,
        software_version="unrecorded experimental worktree",
    ):
        precipitation = load_precipitation_detector(precipitation_detector)

        self.product_type = "conductance_detector"
        self.representation = "detector"
        self.schema_version = SCHEMA_VERSION
        self.conductance_model = CONDUCTANCE_MODEL
        self.conductance_uncertainty_method = (
            CONDUCTANCE_UNCERTAINTY_METHOD
        )
        self.software_version = str(software_version)

        self.precipitation_method = precipitation["method"]
        self.proton_flux_source = precipitation["proton_flux_source"]
        self.proton_energy_model = precipitation["proton_energy_model"]
        self.proton_energy_uncertainty_method = precipitation[
            "proton_energy_uncertainty_method"
        ]
        self.proton_energy_coordinate_note = precipitation[
            "proton_energy_coordinate_note"
        ]
        self.proton_response_energy_min = precipitation[
            "proton_response_energy_min"
        ]
        self.proton_response_energy_max = precipitation[
            "proton_response_energy_max"
        ]
        self.proton_operation_order = precipitation[
            "proton_operation_order"
        ]
        self.count_uncertainty_method = precipitation[
            "count_uncertainty_method"
        ]
        if self.proton_energy_model == "constant":
            self.proton_energy_constant = precipitation[
                "proton_energy_constant"
            ]
            self.proton_energy_uncertainty_constant = precipitation[
                "proton_energy_uncertainty_constant"
            ]

        self.source_precipitation_detector = precipitation["source_file"]
        self.source_precipitation_detector_sha256 = precipitation["source_sha256"]
        self.source_precipitation_detector_size_bytes = precipitation[
            "source_size_bytes"
        ]
        self.source_precipitation_detector_mtime_ns = precipitation[
            "source_mtime_ns"
        ]
        self.source_precipitation_detector_schema_version = precipitation[
            "schema_version"
        ]
        self.source_precipitation_software_version = precipitation[
            "source_software_version"
        ]
        self.source_fuv_detector = precipitation["source_fuv_detector"]
        self.source_fuv_detector_sha256 = precipitation[
            "source_fuv_detector_sha256"
        ]
        self.source_preprocessing_label = precipitation[
            "source_preprocessing_label"
        ]
        self.source_fuv_detector_time_decoding = precipitation[
            "source_fuv_detector_time_decoding"
        ]
        self.coordinate_system = precipitation["coordinate_system"]
        self.reference_height_km = precipitation["reference_height_km"]
        self.kp_provenance = dict(precipitation["kp_provenance"])
        self.shape = precipitation["shape"]

        for name in (
            TIME_FIELDS + TIME_VALUE_FIELDS + GEOMETRY_FIELDS
            + PRECIPITATION_FIELDS
            + (
                "detector_row", "detector_column", "method_quality_weight",
                "method_valid", "Ep_clipping_flag",
            )
        ):
            setattr(self, name, np.asarray(precipitation[name]).copy())

        # Apply the forward model only once on the stored precipitation state.
        result = robinson_conductance(
            self.E0, self.Fe, self.dE0, self.dFe, self.varE0Fe
        )
        for name in ("P", "H", "dP", "dH"):
            values = np.asarray(result[name], dtype=float)
            if values.shape != self.shape:
                raise ValueError(f"{name} must have shape {self.shape}")
            setattr(self, name, values)

        self.conductance_valid = (
            self.method_valid & np.isfinite(self.P) & np.isfinite(self.H)
        )
        self.conductance_uncertainty_valid = (
            self.conductance_valid
            & np.isfinite(self.dP)
            & np.isfinite(self.dH)
        )
        self.P[~self.conductance_valid] = np.nan
        self.H[~self.conductance_valid] = np.nan
        self.dP[~self.conductance_uncertainty_valid] = np.nan
        self.dH[~self.conductance_uncertainty_valid] = np.nan

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
        variable = nc.createVariable(
            name, "f8", ("time",), fill_value=np.nan
        )
        variable[:] = encoded
        variable.units = units
        variable.calendar = "standard"
        variable.time_zone = "UTC"

    def to_nc(self, filename):
        """Write one self-contained detector-space conductance orbit."""

        with Dataset(filename, "w", format="NETCDF4") as nc:
            nc.createDimension("time", self.shape[0])
            nc.createDimension("row", self.shape[1])
            nc.createDimension("column", self.shape[2])

            nc.product_type = self.product_type
            nc.representation = self.representation
            nc.schema_version = self.schema_version
            nc.conductance_model = self.conductance_model
            nc.conductance_uncertainty_method = (
                self.conductance_uncertainty_method
            )
            nc.precipitation_method = self.precipitation_method
            nc.proton_flux_source = self.proton_flux_source
            nc.proton_energy_model = self.proton_energy_model
            nc.proton_energy_uncertainty_method = (
                self.proton_energy_uncertainty_method
            )
            nc.proton_energy_coordinate_note = (
                self.proton_energy_coordinate_note
            )
            nc.proton_response_energy_min = self.proton_response_energy_min
            nc.proton_response_energy_max = self.proton_response_energy_max
            nc.proton_operation_order = self.proton_operation_order
            nc.count_uncertainty_method = self.count_uncertainty_method
            if self.proton_energy_model == "constant":
                nc.proton_energy_constant = self.proton_energy_constant
                nc.proton_energy_uncertainty_constant = (
                    self.proton_energy_uncertainty_constant
                )

            nc.source_precipitation_detector = (
                self.source_precipitation_detector
            )
            nc.source_precipitation_detector_sha256 = (
                self.source_precipitation_detector_sha256
            )
            nc.source_precipitation_detector_size_bytes = np.int64(
                self.source_precipitation_detector_size_bytes
            )
            nc.source_precipitation_detector_mtime_ns = np.int64(
                self.source_precipitation_detector_mtime_ns
            )
            nc.source_precipitation_detector_schema_version = (
                self.source_precipitation_detector_schema_version
            )
            nc.source_precipitation_software_version = (
                self.source_precipitation_software_version
            )
            nc.source_fuv_detector = self.source_fuv_detector
            nc.source_fuv_detector_sha256 = self.source_fuv_detector_sha256
            nc.source_preprocessing_label = self.source_preprocessing_label
            nc.source_fuv_detector_time_decoding = (
                self.source_fuv_detector_time_decoding
            )
            nc.coordinate_system = self.coordinate_system
            nc.reference_height_km = self.reference_height_km
            nc.software_version = self.software_version
            for name, value in self.kp_provenance.items():
                nc.setncattr(f"kp_{name}", value)

            time_units = "seconds since 2000-01-01 00:00:00"
            for name in TIME_FIELDS:
                self._write_time(
                    nc, name, getattr(self, name), time_units
                )

            for name in (
                "wic_source_index", "si12_source_index", "si13_source_index"
            ):
                variable = nc.createVariable(name, "i4", ("time",))
                variable[:] = getattr(self, name)

            for name in (
                "wic_frame_quality", "si12_frame_quality",
                "si13_frame_quality",
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

            for name, units in (("Kp", "1"), ("ssalon", "degrees")):
                variable = nc.createVariable(name, "f4", ("time",))
                variable[:] = getattr(self, name)
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
                "method_quality_weight": (self.method_quality_weight, "1"),
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
                "P": (self.P, "S"),
                "H": (self.H, "S"),
                "dP": (self.dP, "S"),
                "dH": (self.dH, "S"),
            }
            for name, (values, units) in fields.items():
                dtype = "f8" if name in (
                    "glat", "glon", "mlat", "mlon", "mlt"
                ) else "f4"
                variable = nc.createVariable(
                    name, dtype, dimensions, zlib=True
                )
                variable[:] = values
                variable.units = units

            for name in (
                "method_valid", "conductance_valid",
                "conductance_uncertainty_valid", "Ep_clipping_flag",
            ):
                variable = nc.createVariable(
                    name, "i1", dimensions, zlib=True
                )
                variable[:] = getattr(self, name).astype(np.int8)
