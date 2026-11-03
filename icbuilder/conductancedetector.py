"""Calculate conductance on the WIC detector geometry."""

#%% Imports

from functools import lru_cache
from hashlib import sha256
from importlib.metadata import version
from inspect import signature
from pathlib import Path

import numpy as np
from icreader import open_product
from icphysics import (
    FANG_2025_ENERGY_RANGE_KEV,
    FANG_2025_MODEL_DOI,
    fang_2025_conductance_uncertainty,
    robinson_conductance,
)
from netCDF4 import Dataset, date2num

from .fuvdetector import SOURCE_TIME_DECODING, source_identity
from .precipitationdetector import (
    COUNT_UNCERTAINTY_MODE,
    PREPROCESSING_METADATA_ATTRIBUTES,
    PRECIPITATION_METHODS,
    MODELED_ELECTRON_ENERGY_METHODS,
    ELECTRON_METADATA_ATTRIBUTES,
    SCHEMA_VERSION as PRECIPITATION_SCHEMA_VERSION,
)


#%% Product configuration

SCHEMA_VERSION = 3
CONDUCTANCE_MODEL = "robinson"
CONDUCTANCE_MODELS = ("robinson", "fang")
ROBINSON_CONDUCTANCE_UNCERTAINTY_METHOD = (
    "first-order Robinson propagation of Product-2 dE0, dFe, and varE0Fe; "
    "at Fe=0, dP and dH are one-sided excursions from zero to dFe"
)
FANG_CONDUCTANCE_UNCERTAINTY_METHOD = (
    "analytical first-order Fang 2025 propagation of Product-2 dE0, dFe, "
    "and varE0Fe with magnetic-field uncertainty fixed to zero"
)
CONDUCTANCE_UNCERTAINTY_METHOD = ROBINSON_CONDUCTANCE_UNCERTAINTY_METHOD
FANG_LOW_ENERGY_POLICY = "clip"
FANG_MAGNETIC_FIELD_HEIGHT_KM = 150.0
FANG_MAGNETIC_FIELD_CHUNK_SIZE = 65536
FANG_METADATA_ATTRIBUTES = (
    "fang_model_doi",
    "fang_energy_domain_keV",
    "fang_energy_definition",
    "fang_low_energy_policy",
    "fang_magnetic_field_model",
    "fang_magnetic_field_package_version",
    "fang_magnetic_field_coefficients_sha256",
    "fang_magnetic_field_height_km",
    "fang_magnetic_field_epoch_policy",
    "fang_magnetic_field_coordinate_mapping",
    "fang_magnetic_field_uncertainty_method",
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


#%% Fang magnetic field and provenance

def _python_datetime(value):
    """Return a timezone-free Python datetime accepted by Apex and ppigrf."""

    if isinstance(value, np.datetime64):
        return value.astype("datetime64[us]").astype(object)
    if hasattr(value, "to_pydatetime"):
        return value.to_pydatetime()
    return value


def fang_magnetic_field_at_150km(
    time,
    glat,
    glon,
    source_height_km,
    valid=None,
    chunk_size=FANG_MAGNETIC_FIELD_CHUNK_SIZE,
):
    """Evaluate 150 km IGRF above stored geographic detector footprints."""

    import ppigrf

    glat = np.asarray(glat, dtype=float)
    glon = np.asarray(glon, dtype=float)
    if glat.shape != glon.shape or glat.ndim != 3:
        raise ValueError("glat and glon must share (time, row, column) shape")
    if len(time) != glat.shape[0]:
        raise ValueError("time length must match detector frames")
    if not np.isfinite(source_height_km) or source_height_km <= 0:
        raise ValueError("source_height_km must be positive and finite")
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")

    selected = np.isfinite(glat) & np.isfinite(glon)
    if valid is not None:
        valid = np.asarray(valid, dtype=bool)
        if valid.shape != glat.shape:
            raise ValueError("valid must match detector geometry")
        selected &= valid

    latitude = np.where(selected, glat, np.nan)
    longitude = np.where(selected, glon, np.nan)
    field = np.full(glat.shape, np.nan)
    for frame in range(glat.shape[0]):
        indices = np.flatnonzero(selected[frame].ravel())
        if indices.size == 0:
            continue
        epoch = _python_datetime(time[frame])
        frame_latitude = latitude[frame].ravel()
        frame_longitude = longitude[frame].ravel()
        frame_field = field[frame].ravel()

        for start in range(0, indices.size, chunk_size):
            chunk = indices[start:start + chunk_size]
            east, north, up = ppigrf.igrf(
                frame_longitude[chunk],
                frame_latitude[chunk],
                FANG_MAGNETIC_FIELD_HEIGHT_KM,
                epoch,
            )
            magnitude = np.sqrt(
                np.asarray(east).reshape(-1)**2
                + np.asarray(north).reshape(-1)**2
                + np.asarray(up).reshape(-1)**2
            ) / 1000.0
            frame_field[chunk] = magnitude

    return latitude, longitude, field


@lru_cache(maxsize=1)
def fang_conductance_metadata():
    """Return reproducible Fang and IGRF metadata for Product 3."""

    import ppigrf

    coefficients = Path(signature(ppigrf.igrf).parameters["coeff_fn"].default)
    return {
        "fang_model_doi": FANG_2025_MODEL_DOI,
        "fang_energy_domain_keV": (
            f"{FANG_2025_ENERGY_RANGE_KEV[0]:g}--"
            f"{FANG_2025_ENERGY_RANGE_KEV[1]:g}"
        ),
        "fang_energy_definition": (
            "Maxwellian average electron energy; twice Fang characteristic "
            "energy"
        ),
        "fang_low_energy_policy": FANG_LOW_ENERGY_POLICY,
        "fang_magnetic_field_model": coefficients.stem,
        "fang_magnetic_field_package_version": version("ppigrf"),
        "fang_magnetic_field_coefficients_sha256": sha256(
            coefficients.read_bytes()
        ).hexdigest(),
        "fang_magnetic_field_height_km": FANG_MAGNETIC_FIELD_HEIGHT_KM,
        "fang_magnetic_field_epoch_policy": (
            "each IMAGE frame time for both Apex and IGRF"
        ),
        "fang_magnetic_field_coordinate_mapping": (
            "Product-2 stored geographic latitude/longitude at "
            "reference_height_km used directly; IGRF evaluated at 150 km at "
            "the same geographic latitude/longitude"
        ),
        "fang_magnetic_field_uncertainty_method": (
            "not modelled; magnetic_field_std=0"
        ),
    }


#%% Product-2 loading


def load_precipitation_detector(filename):
    """Load the Product-2 state needed by detector conductance."""

    filename = Path(filename)
    with open_product(filename) as source:
        if (
            source.product_type != "precipitation_detector"
            or source.representation != "detector"
            or int(source.schema_version) != PRECIPITATION_SCHEMA_VERSION
            or source.method not in PRECIPITATION_METHODS
            or source.count_uncertainty_mode != COUNT_UNCERTAINTY_MODE
            or source.source_fuv_detector_time_decoding != SOURCE_TIME_DECODING
        ):
            raise ValueError(
                f"{filename} is not a supported schema-"
                f"{PRECIPITATION_SCHEMA_VERSION} supported-method "
                "precipitation_detector product"
            )

        shape = source.shape
        product = {
            "source_file": str(filename),
            "schema_version": int(source.schema_version),
            "method": source.method,
            "proton_flux_source": source.proton_flux_source,
            "proton_energy_model": source.proton_energy_model,
            "proton_energy_uncertainty_method": (
                source.proton_energy_uncertainty_method
            ),
            "proton_energy_coordinate_note": (
                source.proton_energy_coordinate_note
            ),
            "proton_response_energy_min": float(
                source.proton_response_energy_min
            ),
            "proton_response_energy_max": float(
                source.proton_response_energy_max
            ),
            "proton_operation_order": source.proton_operation_order,
            "count_uncertainty_mode": source.count_uncertainty_mode,
            "count_uncertainty_method": source.count_uncertainty_method,
            "coordinate_system": source.coordinate_system,
            "reference_height_km": float(source.reference_height_km),
            "source_software_version": source.software_version,
            "source_fuv_detector": source.source_fuv_detector,
            "source_fuv_detector_sha256": source.source_fuv_detector_sha256,
            "source_preprocessing_label": source.source_preprocessing_label,
            "source_fuv_detector_time_decoding": (
                source.source_fuv_detector_time_decoding
            ),
            "shape": shape,
            "detector_row": np.asarray(source.detector_row, dtype=int).copy(),
            "detector_column": np.asarray(
                source.detector_column, dtype=int
            ).copy(),
            "kp_provenance": dict(source.kp_provenance),
        }
        if source.method in MODELED_ELECTRON_ENERGY_METHODS:
            for name in ELECTRON_METADATA_ATTRIBUTES:
                product[name] = source.attrs[name]
            product["E0_response_range_flag"] = np.asarray(source.read("E0_response_range_flag"), dtype=np.int8)
        for name in PREPROCESSING_METADATA_ATTRIBUTES:
            product[name] = source.attrs[name]
        if product["proton_energy_model"] == "constant":
            product["proton_energy_constant"] = float(source.proton_energy_constant)
            product["proton_energy_uncertainty_constant"] = float(
                source.proton_energy_uncertainty_constant
            )

        for name in TIME_FIELDS:
            product[name] = np.asarray(
                getattr(source, name), dtype=object
            ).copy()

        for name in TIME_VALUE_FIELDS:
            dtype = int if name.endswith("_index") or name.endswith(
                "_frame_quality"
            ) else float
            product[name] = np.asarray(
                getattr(source, name), dtype=dtype
            ).copy()

        for name in GEOMETRY_FIELDS + PRECIPITATION_FIELDS + (
            "method_quality_weight",
        ):
            product[name] = np.asarray(source.read(name), dtype=float)

        for name in ("method_valid", "Ep_clipping_flag"):
            product[name] = np.asarray(source.read(name), dtype=bool)

    product.update(source_identity(filename))
    return product


#%% Detector conductance

class ConductanceDetector:
    """Apply a selected conductance model to one detector Product-2 orbit."""

    def __init__(
        self,
        precipitation_detector,
        *,
        conductance_model=CONDUCTANCE_MODEL,
        software_version="unrecorded experimental worktree",
    ):
        if conductance_model not in CONDUCTANCE_MODELS:
            raise ValueError(
                f"conductance_model must be one of {CONDUCTANCE_MODELS}"
            )
        precipitation = load_precipitation_detector(precipitation_detector)

        self.product_type = "conductance_detector"
        self.representation = "detector"
        self.schema_version = SCHEMA_VERSION
        self.conductance_model = conductance_model
        self.conductance_uncertainty_method = (
            ROBINSON_CONDUCTANCE_UNCERTAINTY_METHOD
            if conductance_model == "robinson"
            else FANG_CONDUCTANCE_UNCERTAINTY_METHOD
        )
        self.software_version = str(software_version)

        self.precipitation_method = precipitation["method"]
        if self.precipitation_method in MODELED_ELECTRON_ENERGY_METHODS:
            for name in ELECTRON_METADATA_ATTRIBUTES:
                setattr(self, name, precipitation[name])
            self.E0_response_range_flag = precipitation["E0_response_range_flag"].copy()
        self.count_source = str(precipitation["count_source"])
        self.smoothed = bool(int(precipitation["smoothed"]))
        self.smoothing_method = str(precipitation["smoothing_method"])
        for sensor in ("wic", "si12", "si13"):
            setattr(
                self,
                f"{sensor}_smoothing_width_pixels",
                float(precipitation[f"{sensor}_smoothing_width_pixels"]),
            )
            setattr(
                self,
                f"{sensor}_smoothing_applied",
                bool(int(precipitation[f"{sensor}_smoothing_applied"])),
            )
        self.smoothing_width_definition = str(
            precipitation["smoothing_width_definition"]
        )
        self.smoothing_operation_order = str(
            precipitation["smoothing_operation_order"]
        )
        self.smoothing_variance_method = str(
            precipitation["smoothing_variance_method"]
        )
        self.method_quality_weight_method = str(
            precipitation["method_quality_weight_method"]
        )
        self.method_quality_weight_floor = float(
            precipitation["method_quality_weight_floor"]
        )
        self.method_quality_weight_spatial_propagation = str(
            precipitation["method_quality_weight_spatial_propagation"]
        )
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
        self.count_uncertainty_mode = precipitation[
            "count_uncertainty_mode"
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

        # Apply the selected forward model once on the stored precipitation.
        if conductance_model == "robinson":
            result = robinson_conductance(
                self.E0, self.Fe, self.dE0, self.dFe, self.varE0Fe
            )
            output = {
                "P": result["P"],
                "H": result["H"],
                "dP": result["dP"],
                "dH": result["dH"],
            }
            uncertainty_valid = (
                np.isfinite(output["dP"]) & np.isfinite(output["dH"])
            )
        else:
            for name, value in fang_conductance_metadata().items():
                setattr(self, name, value)
            (
                self.fang_geographic_latitude_at_reference_height,
                self.fang_geographic_longitude_at_reference_height,
                self.fang_magnetic_field_150km_microtesla,
            ) = fang_magnetic_field_at_150km(
                self.time,
                self.glat,
                self.glon,
                self.reference_height_km,
                valid=self.method_valid,
            )
            result = fang_2025_conductance_uncertainty(
                self.E0,
                self.Fe,
                self.fang_magnetic_field_150km_microtesla,
                mean_energy_std=self.dE0,
                energy_flux_std=self.dFe,
                magnetic_field_std=np.zeros(self.shape),
                energy_flux_covariance=self.varE0Fe,
                low_energy_policy=FANG_LOW_ENERGY_POLICY,
            )
            output = {
                "P": result["P"],
                "H": result["H"],
                "dP": result["P_std"],
                "dH": result["H_std"],
            }
            uncertainty_valid = result["uncertainty_valid"]
            self.fang_evaluated_energy = np.asarray(
                result["evaluated_energy"], dtype=float
            )
            self.fang_within_model_domain = np.asarray(
                result["within_model_domain"], dtype=bool
            )
            self.fang_energy_clipped_low = np.asarray(
                result["energy_clipped_low"], dtype=bool
            )
            self.fang_uncertainty_straddles_low_energy_boundary = np.asarray(
                result["uncertainty_straddles_low_energy_boundary"],
                dtype=bool,
            )

        for name in ("P", "H", "dP", "dH"):
            values = np.asarray(output[name], dtype=float)
            if values.shape != self.shape:
                raise ValueError(f"{name} must have shape {self.shape}")
            setattr(self, name, values)

        self.conductance_valid = (
            self.method_valid & np.isfinite(self.P) & np.isfinite(self.H)
        )
        self.conductance_uncertainty_valid = (
            self.conductance_valid
            & uncertainty_valid
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
            if self.conductance_model == "fang":
                for name in FANG_METADATA_ATTRIBUTES:
                    nc.setncattr(name, getattr(self, name))
            nc.precipitation_method = self.precipitation_method
            if self.precipitation_method in MODELED_ELECTRON_ENERGY_METHODS:
                for name in ELECTRON_METADATA_ATTRIBUTES:
                    nc.setncattr(name, getattr(self, name))
            nc.count_source = self.count_source
            nc.smoothed = np.int8(self.smoothed)
            nc.smoothing_method = self.smoothing_method
            nc.wic_smoothing_width_pixels = self.wic_smoothing_width_pixels
            nc.si12_smoothing_width_pixels = self.si12_smoothing_width_pixels
            nc.si13_smoothing_width_pixels = self.si13_smoothing_width_pixels
            for sensor in ("wic", "si12", "si13"):
                nc.setncattr(
                    f"{sensor}_smoothing_applied",
                    np.int8(getattr(self, f"{sensor}_smoothing_applied")),
                )
            nc.smoothing_width_definition = self.smoothing_width_definition
            nc.smoothing_operation_order = self.smoothing_operation_order
            nc.smoothing_variance_method = self.smoothing_variance_method
            nc.method_quality_weight_method = self.method_quality_weight_method
            nc.method_quality_weight_floor = self.method_quality_weight_floor
            nc.method_quality_weight_spatial_propagation = (
                self.method_quality_weight_spatial_propagation
            )
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
            nc.count_uncertainty_mode = self.count_uncertainty_mode
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
            if self.conductance_model == "fang":
                fields.update({
                    "fang_geographic_latitude_at_reference_height": (
                        self.fang_geographic_latitude_at_reference_height,
                        "degrees_north",
                    ),
                    "fang_geographic_longitude_at_reference_height": (
                        self.fang_geographic_longitude_at_reference_height,
                        "degrees_east",
                    ),
                    "fang_magnetic_field_150km_microtesla": (
                        self.fang_magnetic_field_150km_microtesla,
                        "microtesla",
                    ),
                    "fang_evaluated_energy": (
                        self.fang_evaluated_energy,
                        "keV",
                    ),
                })
            for name, (values, units) in fields.items():
                dtype = "f8" if name in (
                    "glat", "glon", "mlat", "mlon", "mlt"
                ) else "f4"
                variable = nc.createVariable(
                    name, dtype, dimensions, zlib=True
                )
                variable[:] = values
                variable.units = units

            if self.precipitation_method in MODELED_ELECTRON_ENERGY_METHODS:
                flag = nc.createVariable("E0_response_range_flag", "i1", dimensions, zlib=True)
                flag[:] = self.E0_response_range_flag
                flag.flag_values = np.asarray([-1, 0, 1, 2], dtype=np.int8)
                flag.flag_meanings = "unavailable in_range below above"

            for name in (
                "method_valid", "conductance_valid",
                "conductance_uncertainty_valid", "Ep_clipping_flag",
            ):
                variable = nc.createVariable(
                    name, "i1", dimensions, zlib=True
                )
                variable[:] = getattr(self, name).astype(np.int8)

            if self.conductance_model == "fang":
                fang_flags = {
                    "fang_within_model_domain": (
                        "original energy is within the Fang 1--20 keV domain"
                    ),
                    "fang_energy_clipped_low": (
                        "original energy is positive and below 1 keV"
                    ),
                    "fang_uncertainty_straddles_low_energy_boundary": (
                        "one-standard-deviation energy interval crosses 1 keV"
                    ),
                }
                for name, description in fang_flags.items():
                    variable = nc.createVariable(
                        name, "i1", dimensions, zlib=True
                    )
                    variable[:] = getattr(self, name).astype(np.int8)
                    variable.long_name = description
