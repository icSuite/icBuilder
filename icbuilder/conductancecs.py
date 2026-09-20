"""Reduce detector-space conductance onto the frozen Cubed-Sphere grid."""

#%% Imports

from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from .conductancedetector import (
    CONDUCTANCE_MODEL,
    SCHEMA_VERSION as DETECTOR_SCHEMA_VERSION,
)
from .detectorcs import (
    BINNING_METHOD,
    UNCERTAINTY_METHOD,
    read_common_time_fields,
    read_variable,
    reduce_area_mean,
    reduce_covariance,
    reduce_flag_fraction,
    reduce_measurement_variance,
    reduce_support,
    write_common_cs_coordinates,
)
from .fuvdetector import source_identity
from .grids import DETECTOR_CS_COORDINATE_SHA256, DETECTOR_CS_GRID_ID
from .precipitationdetector import COUNT_UNCERTAINTY_MODE


#%% Product contract

SCHEMA_VERSION = 1

MEAN_FIELDS = {
    "method_quality_weight": "1",
    "E0": "keV",
    "Fe": "mW m-2",
    "P": "S",
    "H": "S",
}
UNCERTAINTY_FIELDS = {
    "dE0": "keV",
    "dFe": "mW m-2",
    "dP": "S",
    "dH": "S",
}
COVARIANCE_FIELDS = {
    "varE0Fe": "keV mW m-2",
}
RESULT_FIELDS = tuple(
    list(MEAN_FIELDS)
    + list(UNCERTAINTY_FIELDS)
    + list(COVARIANCE_FIELDS)
    + [
        "source_count", "coverage",
        "uncertainty_source_count", "uncertainty_coverage",
        "conductance_valid", "conductance_uncertainty_valid",
        "Ep_clipping_fraction", "Ep_clipping_any",
    ]
)


def _require_source(source, filename):
    """Validate the detector Product-3 boundary used by this representation."""

    if (
        source.product_type != "conductance_detector"
        or source.representation != "detector"
        or int(source.schema_version) != DETECTOR_SCHEMA_VERSION
        or source.conductance_model != CONDUCTANCE_MODEL
        or source.count_uncertainty_mode != COUNT_UNCERTAINTY_MODE
    ):
        raise ValueError(
            f"{filename} is not a supported schema-{DETECTOR_SCHEMA_VERSION} "
            f"{CONDUCTANCE_MODEL} conductance_detector product"
        )

    shape = (
        len(source.dimensions["time"]),
        len(source.dimensions["row"]),
        len(source.dimensions["column"]),
    )
    if any(length == 0 for length in shape):
        raise ValueError("conductance_detector dimensions must be non-empty")
    for name in (
        "mlat", "mlt", "conductance_valid",
        "conductance_uncertainty_valid", "Ep_clipping_flag",
        *MEAN_FIELDS, *UNCERTAINTY_FIELDS, *COVARIANCE_FIELDS,
    ):
        if source.variables[name].shape != shape:
            raise ValueError(f"{name} must have shape {shape}")
    return shape


#%% Product implementation

class ConductanceCS:
    """Area-reduced representation of one detector Product-3 orbit."""

    def __init__(self, filename, *, grid, software_version):
        filename = Path(filename)
        self.source_identity = source_identity(filename)
        self.grid = grid
        self.grid_id = DETECTOR_CS_GRID_ID
        self.grid_coordinate_sha256 = DETECTOR_CS_COORDINATE_SHA256
        self.software_version = str(software_version)
        self.companion_precipitation_cs = ""

        with Dataset(filename) as source:
            detector_shape = _require_source(source, filename)
            self.shape = (detector_shape[0], *grid.shape)
            self.reference_height_km = float(source.reference_height_km)
            self.time_fields = read_common_time_fields(source)

            self.conductance_model = source.conductance_model
            self.conductance_uncertainty_method = (
                source.conductance_uncertainty_method
            )
            self.precipitation_method = source.precipitation_method
            self.proton_flux_source = source.proton_flux_source
            self.proton_energy_model = source.proton_energy_model
            self.proton_energy_uncertainty_method = (
                source.proton_energy_uncertainty_method
            )
            self.proton_energy_coordinate_note = (
                source.proton_energy_coordinate_note
            )
            self.proton_response_energy_min = float(
                source.proton_response_energy_min
            )
            self.proton_response_energy_max = float(
                source.proton_response_energy_max
            )
            self.proton_operation_order = source.proton_operation_order
            self.count_uncertainty_mode = source.count_uncertainty_mode
            self.count_uncertainty_method = source.count_uncertainty_method
            self.coordinate_system = source.coordinate_system
            self.source_software_version = source.software_version
            self.source_precipitation_detector = (
                source.source_precipitation_detector
            )
            self.source_precipitation_detector_sha256 = (
                source.source_precipitation_detector_sha256
            )
            self.source_fuv_detector = source.source_fuv_detector
            self.source_fuv_detector_sha256 = source.source_fuv_detector_sha256
            self.source_preprocessing_label = source.source_preprocessing_label
            self.source_fuv_detector_time_decoding = (
                source.source_fuv_detector_time_decoding
            )
            self.kp_provenance = {
                name[3:]: source.getncattr(name)
                for name in source.ncattrs()
                if name.startswith("kp_")
            }
            if self.proton_energy_model == "constant":
                self.proton_energy_constant = float(
                    source.proton_energy_constant
                )
                self.proton_energy_uncertainty_constant = float(
                    source.proton_energy_uncertainty_constant
                )

        for name in MEAN_FIELDS | UNCERTAINTY_FIELDS | COVARIANCE_FIELDS:
            setattr(self, name, np.full(self.shape, np.nan))
        self.source_count = np.zeros(self.shape, dtype=np.int32)
        self.coverage = np.zeros(self.shape)
        self.uncertainty_source_count = np.zeros(self.shape, dtype=np.int32)
        self.uncertainty_coverage = np.zeros(self.shape)
        self.conductance_valid = np.zeros(self.shape, dtype=bool)
        self.conductance_uncertainty_valid = np.zeros(self.shape, dtype=bool)
        self.Ep_clipping_fraction = np.full(self.shape, np.nan)
        self.Ep_clipping_any = np.zeros(self.shape, dtype=bool)

    def reduce_frame(self, source, frame, mapping, cell_area):
        """Reduce one detector frame with a caller-supplied common mapping."""

        output_shape = self.grid.shape
        conductance_valid = read_variable(
            source.variables["conductance_valid"], frame, bool
        )

        for name in MEAN_FIELDS:
            values = read_variable(source.variables[name], frame)
            reduced, _ = reduce_area_mean(
                values, conductance_valid, mapping, output_shape
            )
            getattr(self, name)[frame] = reduced

        for name in UNCERTAINTY_FIELDS:
            uncertainty = read_variable(source.variables[name], frame)
            variance, _ = reduce_measurement_variance(
                uncertainty**2, conductance_valid, mapping, output_shape
            )
            getattr(self, name)[frame] = np.sqrt(variance)

        for name in COVARIANCE_FIELDS:
            covariance = read_variable(source.variables[name], frame)
            reduced, _ = reduce_covariance(
                covariance, conductance_valid, mapping, output_shape
            )
            getattr(self, name)[frame] = reduced

        count, coverage = reduce_support(
            conductance_valid, mapping, cell_area, output_shape
        )
        uncertainty_valid = read_variable(
            source.variables["conductance_uncertainty_valid"], frame, bool
        )
        uncertainty_count, uncertainty_coverage = reduce_support(
            uncertainty_valid, mapping, cell_area, output_shape
        )
        clipping_fraction, clipping_any = reduce_flag_fraction(
            read_variable(source.variables["Ep_clipping_flag"], frame, bool),
            conductance_valid,
            mapping,
            output_shape,
        )

        self.source_count[frame] = count
        self.coverage[frame] = coverage
        self.uncertainty_source_count[frame] = uncertainty_count
        self.uncertainty_coverage[frame] = uncertainty_coverage
        self.conductance_valid[frame] = coverage > 0
        self.conductance_uncertainty_valid[frame] = uncertainty_coverage > 0
        self.Ep_clipping_fraction[frame] = clipping_fraction
        self.Ep_clipping_any[frame] = clipping_any

    def reduce(self, source, mappings, cell_area):
        """Reduce the orbit while reading each compressed variable only once."""

        output_shape = self.grid.shape
        conductance_valid = read_variable(
            source.variables["conductance_valid"], dtype=bool
        )
        uncertainty_valid = read_variable(
            source.variables["conductance_uncertainty_valid"], dtype=bool
        )

        for name in MEAN_FIELDS:
            values = read_variable(source.variables[name])
            for frame, mapping in enumerate(mappings):
                reduced, _ = reduce_area_mean(
                    values[frame],
                    conductance_valid[frame],
                    mapping,
                    output_shape,
                )
                getattr(self, name)[frame] = reduced

        for name in UNCERTAINTY_FIELDS:
            uncertainty = read_variable(source.variables[name])
            for frame, mapping in enumerate(mappings):
                variance, _ = reduce_measurement_variance(
                    uncertainty[frame] ** 2,
                    conductance_valid[frame],
                    mapping,
                    output_shape,
                )
                getattr(self, name)[frame] = np.sqrt(variance)

        for name in COVARIANCE_FIELDS:
            covariance = read_variable(source.variables[name])
            for frame, mapping in enumerate(mappings):
                reduced, _ = reduce_covariance(
                    covariance[frame],
                    conductance_valid[frame],
                    mapping,
                    output_shape,
                )
                getattr(self, name)[frame] = reduced

        clipping_flag = read_variable(
            source.variables["Ep_clipping_flag"], dtype=bool
        )
        for frame, mapping in enumerate(mappings):
            count, coverage = reduce_support(
                conductance_valid[frame], mapping, cell_area, output_shape
            )
            uncertainty_count, uncertainty_coverage = reduce_support(
                uncertainty_valid[frame], mapping, cell_area, output_shape
            )
            clipping_fraction, clipping_any = reduce_flag_fraction(
                clipping_flag[frame],
                conductance_valid[frame],
                mapping,
                output_shape,
            )
            self.source_count[frame] = count
            self.coverage[frame] = coverage
            self.uncertainty_source_count[frame] = uncertainty_count
            self.uncertainty_coverage[frame] = uncertainty_coverage
            self.conductance_valid[frame] = coverage > 0
            self.conductance_uncertainty_valid[frame] = (
                uncertainty_coverage > 0
            )
            self.Ep_clipping_fraction[frame] = clipping_fraction
            self.Ep_clipping_any[frame] = clipping_any

    def to_nc(self, filename):
        """Write one fixed-grid conductance orbit."""

        with Dataset(filename, "w", format="NETCDF4") as nc:
            write_common_cs_coordinates(nc, self)

            nc.product_type = "conductance_cs"
            nc.representation = "cs"
            nc.schema_version = SCHEMA_VERSION
            nc.grid_id = self.grid_id
            nc.grid_coordinate_sha256 = self.grid_coordinate_sha256
            nc.binning_method = BINNING_METHOD
            nc.uncertainty_method = UNCERTAINTY_METHOD
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
            nc.proton_energy_coordinate_note = self.proton_energy_coordinate_note
            nc.proton_response_energy_min = self.proton_response_energy_min
            nc.proton_response_energy_max = self.proton_response_energy_max
            nc.proton_operation_order = self.proton_operation_order
            nc.count_uncertainty_mode = self.count_uncertainty_mode
            nc.count_uncertainty_method = self.count_uncertainty_method
            nc.coordinate_system = self.coordinate_system
            nc.reference_height_km = self.reference_height_km
            nc.software_version = self.software_version
            nc.source_conductance_detector = self.source_identity["source_file"]
            nc.source_conductance_detector_sha256 = self.source_identity[
                "source_sha256"
            ]
            nc.source_conductance_detector_size_bytes = np.int64(
                self.source_identity["source_size_bytes"]
            )
            nc.source_conductance_detector_mtime_ns = np.int64(
                self.source_identity["source_mtime_ns"]
            )
            nc.source_conductance_detector_schema_version = DETECTOR_SCHEMA_VERSION
            nc.source_conductance_software_version = self.source_software_version
            nc.source_precipitation_detector = (
                self.source_precipitation_detector
            )
            nc.source_precipitation_detector_sha256 = (
                self.source_precipitation_detector_sha256
            )
            nc.companion_precipitation_cs = self.companion_precipitation_cs
            nc.source_fuv_detector = self.source_fuv_detector
            nc.source_fuv_detector_sha256 = self.source_fuv_detector_sha256
            nc.source_preprocessing_label = self.source_preprocessing_label
            nc.source_fuv_detector_time_decoding = (
                self.source_fuv_detector_time_decoding
            )
            if self.proton_energy_model == "constant":
                nc.proton_energy_constant = self.proton_energy_constant
                nc.proton_energy_uncertainty_constant = (
                    self.proton_energy_uncertainty_constant
                )
            for name, value in self.kp_provenance.items():
                nc.setncattr(f"kp_{name}", value)

            dimensions = ("time", "dim1", "dim2")
            for name, units in (
                list(MEAN_FIELDS.items())
                + list(UNCERTAINTY_FIELDS.items())
                + list(COVARIANCE_FIELDS.items())
            ):
                variable = nc.createVariable(
                    name, "f4", dimensions, zlib=True
                )
                variable[:] = getattr(self, name)
                variable.units = units

            for name in ("source_count", "uncertainty_source_count"):
                variable = nc.createVariable(name, "i4", dimensions, zlib=True)
                variable[:] = getattr(self, name)
                variable.units = "1"

            for name in (
                "coverage", "uncertainty_coverage", "Ep_clipping_fraction"
            ):
                variable = nc.createVariable(name, "f4", dimensions, zlib=True)
                variable[:] = getattr(self, name)
                variable.units = "1"

            for name in (
                "conductance_valid", "conductance_uncertainty_valid",
                "Ep_clipping_any",
            ):
                variable = nc.createVariable(name, "i1", dimensions, zlib=True)
                variable[:] = getattr(self, name).astype(np.int8)

            nc.coverage_definition = (
                "Valid detector-footprint overlap area divided by projected "
                "target-cell area, clipped to one."
            )
            nc.source_count_definition = (
                "Number of conductance-valid WIC detector footprints "
                "intersecting the cell."
            )
            nc.nonlinear_ordering = (
                "P and H are overlap means of detector-space conductance; they "
                "are not recalculated from binned E0 and Fe."
            )
