"""Calculate precipitation on the WIC detector geometry."""

#%% Imports

from pathlib import Path

import numpy as np
from icreader import open_product
from icphysics import (
    PROTON_RESPONSE_ENERGY_RANGE,
    hardy_ion_precipitation,
    precipitation_from_ratio,
    proton_correct_images,
)
from netCDF4 import Dataset, date2num
from scipy.ndimage import convolve

from .kp import load_gfz_kp, match_gfz_kp
from .fuvdetector import (
    PREPROCESSING_LABEL,
    SCHEMA_VERSION as FUV_SCHEMA_VERSION,
    SOURCE_TIME_DECODING,
    source_identity,
)


#%% Product configuration

SCHEMA_VERSION = 3
PRECIPITATION_METHOD = "image_ratio"
PROTON_ENERGY_MODELS = ("hardy", "constant")
COUNT_SOURCES = ("background_subtracted", "unsubtracted")
SPATIAL_SMOOTHING_KERNELS = ("none", "gaussian", "boxcar")
PROTON_FLUX_SOURCE = "SI12"
COUNT_UNCERTAINTY_MODE = "measurement"
PROTON_OPERATION_ORDER = (
    "SI12 counts are coregistered onto WIC pixels in fuv_detector, then "
    "converted to proton flux on the WIC detector geometry"
)
COUNT_UNCERTAINTY_METHOD = (
    "square roots of Product-1 full detector measurement variances are "
    "supplied to icPhysics without adding a second raw-count Poisson term. "
    "SI variances already include independent-pixel propagation through "
    "coregistration; covariance and background-model uncertainty are not "
    "included"
)
HARDY_COORDINATE_NOTE = (
    "Hardy corrected geomagnetic coordinates approximated by Product-1 "
    "Modified Apex latitude and MLT at 130 km"
)
SPATIAL_SMOOTHING_ORDER = (
    "WIC and SI13 are smoothed independently on their Product-1 valid "
    "support after count-source selection and before SI12 proton correction"
)
SPATIAL_SMOOTHING_VARIANCE_METHOD = (
    "independent-pixel input variances propagated with squared normalized "
    "spatial-kernel weights; smoothing-induced output covariance is not stored"
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
UNSUBTRACTED_FRAME_FIELDS = (
    "wic_unsubtracted_counts", "si12_unsubtracted_counts",
    "si13_unsubtracted_counts", "wic_unsubtracted_valid",
    "si12_unsubtracted_valid", "si13_unsubtracted_valid",
)
TIME_FIELDS = (
    "time", "wic_source_time", "si12_source_time", "si13_source_time",
    "wic_source_index", "si12_source_index", "si13_source_index", "ssalon",
    "wic_frame_quality", "si12_frame_quality", "si13_frame_quality",
)


#%% Product-1 loading


def load_fuv_detector(filename, include_unsubtracted=False):
    """Load the Product-1 fields needed by detector precipitation."""

    filename = Path(filename)
    with open_product(filename) as source:
        if (
            source.product_type != "fuv_detector"
            or source.representation != "detector"
            or int(source.schema_version) != FUV_SCHEMA_VERSION
            or source.preprocessing_label != PREPROCESSING_LABEL
            or source.source_time_decoding != SOURCE_TIME_DECODING
        ):
            raise ValueError(
                f"{filename} is not a supported schema-{FUV_SCHEMA_VERSION} "
                f"{PREPROCESSING_LABEL} fuv_detector product"
            )

        product = {
            "source_file": str(filename),
            "schema_version": int(source.schema_version),
            "preprocessing_label": source.preprocessing_label,
            "source_time_decoding": source.source_time_decoding,
            "source_software_version": source.software_version,
            "coordinate_system": source.coordinate_system,
            "reference_height_km": float(source.reference_height_km),
            "shape": source.shape,
            "detector_row": np.asarray(source.detector_row, dtype=int).copy(),
            "detector_column": np.asarray(
                source.detector_column, dtype=int
            ).copy(),
        }
        for name in TIME_FIELDS:
            if name.endswith("_time") or name == "time":
                product[name] = np.asarray(
                    getattr(source, name), dtype=object
                ).copy()
            elif name.endswith("_index") or name.endswith("_frame_quality"):
                product[name] = np.asarray(
                    getattr(source, name), dtype=int
                ).copy()
            else:
                product[name] = np.asarray(
                    getattr(source, name), dtype=float
                ).copy()
        for name in FRAME_FIELDS:
            dtype = bool if name.endswith("_valid") else float
            product[name] = np.asarray(source.read(name), dtype=dtype)
        if include_unsubtracted:
            if int(source.attrs.get("unsubtracted_counts_stored", 0)) != 1:
                raise ValueError(
                    f"{filename} does not contain unsubtracted Product-1 counts"
                )
            for name in UNSUBTRACTED_FRAME_FIELDS:
                dtype = bool if name.endswith("_valid") else float
                product[name] = np.asarray(source.read(name), dtype=dtype)

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


#%% Diagnostic spatial smoothing

def resolve_smoothing_configuration(kernel, wic_width=None, si13_width=None):
    """Validate smoothing settings and supply method-specific defaults."""

    if kernel not in SPATIAL_SMOOTHING_KERNELS:
        raise ValueError(
            f"spatial_smoothing_kernel must be one of {SPATIAL_SMOOTHING_KERNELS}"
        )

    if kernel == "none":
        supplied = [width for width in (wic_width, si13_width) if width is not None]
        if any(float(width) != 0 for width in supplied):
            raise ValueError("smoothing widths require gaussian or boxcar smoothing")
        return 0.0, 0.0

    default = 1.0 if kernel == "gaussian" else 3.0
    widths = [default if width is None else float(width) for width in (wic_width, si13_width)]
    if any(not np.isfinite(width) or width <= 0 for width in widths):
        raise ValueError("smoothing widths must be positive and finite")
    if kernel == "boxcar" and any(
        not width.is_integer() or int(width) % 2 != 1 for width in widths
    ):
        raise ValueError("boxcar smoothing widths must be positive odd integers")
    return tuple(widths)


def spatial_smoothing_kernel(kernel, width):
    """Return one normalized two-dimensional diagnostic smoothing kernel."""

    if kernel == "gaussian":
        radius = max(1, int(np.ceil(4 * width)))
        coordinate = np.arange(-radius, radius + 1, dtype=float)
        row, column = np.meshgrid(coordinate, coordinate, indexing="ij")
        weights = np.exp(-(row**2 + column**2) / (2 * width**2))
    elif kernel == "boxcar":
        weights = np.ones((int(width), int(width)), dtype=float)
    else:
        raise ValueError("a spatial kernel is only defined for gaussian or boxcar")
    return weights / np.sum(weights)


def smooth_detector_counts(counts, variance, valid, kernel, width):
    """Smooth detector frames and propagate independent-pixel variances."""

    counts = np.asarray(counts, dtype=float)
    variance = np.asarray(variance, dtype=float)
    valid = np.asarray(valid, dtype=bool)
    if counts.shape != variance.shape or counts.shape != valid.shape:
        raise ValueError("counts, variance, and valid must have the same shape")
    if counts.ndim != 3:
        raise ValueError("detector smoothing expects time, row, column arrays")

    weights = spatial_smoothing_kernel(kernel, width)
    smoothed = np.full(counts.shape, np.nan)
    smoothed_variance = np.full(counts.shape, np.nan)
    smoothed_valid = np.zeros(counts.shape, dtype=bool)

    for frame in range(counts.shape[0]):
        count_valid = valid[frame] & np.isfinite(counts[frame])
        normalization = convolve(
            count_valid.astype(float), weights, mode="constant", cval=0.0
        )
        numerator = convolve(
            np.where(count_valid, counts[frame], 0.0),
            weights,
            mode="constant",
            cval=0.0,
        )
        output_valid = count_valid & (normalization > 0)
        smoothed[frame] = np.divide(
            numerator,
            normalization,
            out=np.full(counts.shape[1:], np.nan),
            where=output_valid,
        )

        variance_valid = count_valid & np.isfinite(variance[frame])
        missing_variance_weight = convolve(
            (count_valid & ~variance_valid).astype(float),
            weights,
            mode="constant",
            cval=0.0,
        )
        variance_numerator = convolve(
            np.where(variance_valid, variance[frame], 0.0),
            weights**2,
            mode="constant",
            cval=0.0,
        )
        uncertainty_valid = output_valid & (missing_variance_weight < 1e-12)
        smoothed_variance[frame] = np.divide(
            variance_numerator,
            normalization**2,
            out=np.full(counts.shape[1:], np.nan),
            where=uncertainty_valid,
        )
        smoothed_valid[frame] = output_valid

    return smoothed, smoothed_variance, smoothed_valid


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
        count_source="background_subtracted",
        spatial_smoothing_kernel="none",
        wic_smoothing_width=None,
        si13_smoothing_width=None,
        software_version="unrecorded experimental worktree",
    ):
        if count_source not in COUNT_SOURCES:
            raise ValueError(f"count_source must be one of {COUNT_SOURCES}")
        fuv = load_fuv_detector(
            fuv_detector,
            include_unsubtracted=count_source == "unsubtracted",
        )
        if kp_series is None:
            kp_series = load_gfz_kp()

        self.product_type = "precipitation_detector"
        self.representation = "detector"
        self.schema_version = SCHEMA_VERSION
        self.method = PRECIPITATION_METHOD
        self.count_source = count_source
        self.proton_flux_source = PROTON_FLUX_SOURCE
        self.proton_energy_model = proton_energy_model
        self.proton_energy_constant = float(proton_energy)
        self.proton_energy_uncertainty_constant = float(
            proton_energy_uncertainty
        )
        self.proton_operation_order = PROTON_OPERATION_ORDER
        self.count_uncertainty_mode = COUNT_UNCERTAINTY_MODE
        self.count_uncertainty_method = COUNT_UNCERTAINTY_METHOD
        smoothing_widths = resolve_smoothing_configuration(
            spatial_smoothing_kernel,
            wic_smoothing_width,
            si13_smoothing_width,
        )
        self.spatial_smoothing_kernel = spatial_smoothing_kernel
        self.wic_smoothing_width_pixels = smoothing_widths[0]
        self.si13_smoothing_width_pixels = smoothing_widths[1]
        if self.spatial_smoothing_kernel == "gaussian":
            self.spatial_smoothing_width_definition = (
                "Gaussian sigma in native WIC detector pixels"
            )
        elif self.spatial_smoothing_kernel == "boxcar":
            self.spatial_smoothing_width_definition = (
                "odd square side length in native WIC detector pixels"
            )
        else:
            self.spatial_smoothing_width_definition = "not applicable"
        self.spatial_smoothing_operation_order = SPATIAL_SMOOTHING_ORDER
        self.spatial_smoothing_variance_method = (
            SPATIAL_SMOOTHING_VARIANCE_METHOD
        )
        if self.spatial_smoothing_kernel != "none":
            self.count_uncertainty_method += (
                "; " + self.spatial_smoothing_variance_method
            )
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

        if self.count_source == "unsubtracted":
            for name in UNSUBTRACTED_FRAME_FIELDS:
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

        # 3. Select one count stage on common Product-1 detector support.
        if self.count_source == "background_subtracted":
            wic_counts = self.wic_counts
            si12_counts = self.si12_counts
            si13_counts = self.si13_counts
            self.method_quality_weight_method = (
                "product of Product-1 fuvpy dgweight fields"
            )
        else:
            wic_counts = self.wic_unsubtracted_counts
            si12_counts = self.si12_unsubtracted_counts
            si13_counts = self.si13_unsubtracted_counts
            self.wic_valid &= self.wic_unsubtracted_valid
            self.si12_valid &= self.si12_unsubtracted_valid
            self.si13_valid &= self.si13_unsubtracted_valid
            self.method_quality_weight_method = (
                "uniform weight on successful common detector support; "
                "background-fit dgweight is not used"
            )

        wic_method_valid = self.wic_valid
        si13_method_valid = self.si13_valid
        wic_variance = self.wic_variance
        si13_variance = self.si13_variance

        if self.spatial_smoothing_kernel != "none":
            smoothed = smooth_detector_counts(
                wic_counts,
                self.wic_variance,
                self.wic_valid,
                self.spatial_smoothing_kernel,
                self.wic_smoothing_width_pixels,
            )
            self.wic_smoothed, wic_variance, wic_method_valid = smoothed
            smoothed = smooth_detector_counts(
                si13_counts,
                self.si13_variance,
                self.si13_valid,
                self.spatial_smoothing_kernel,
                self.si13_smoothing_width_pixels,
            )
            self.si13_smoothed, si13_variance, si13_method_valid = smoothed
            self.dwic_smoothed = np.sqrt(wic_variance)
            self.dsi13_smoothed = np.sqrt(si13_variance)
            wic_counts = self.wic_smoothed
            si13_counts = self.si13_smoothed

        input_valid = (
            wic_method_valid
            & self.si12_valid
            & si13_method_valid
            & np.isfinite(self.Ep)
        )
        wic = np.where(input_valid, wic_counts, np.nan)
        si12 = np.where(input_valid, si12_counts, np.nan)
        si13 = np.where(input_valid, si13_counts, np.nan)

        # 4. Infer proton flux from mapped SI12, then correct WIC and SI13.
        dwic = np.where(input_valid, np.sqrt(wic_variance), np.nan)
        dsi12 = np.where(input_valid, np.sqrt(self.si12_variance), np.nan)
        dsi13 = np.where(input_valid, np.sqrt(si13_variance), np.nan)
        self.si12 = si12
        self.dsi12 = dsi12
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
                uncertainty_mode=self.count_uncertainty_mode,
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
        if self.count_source == "background_subtracted":
            method_quality_weight = (
                self.wic_quality_weight
                * self.si12_quality_weight
                * self.si13_quality_weight
            )
        else:
            method_quality_weight = np.ones(self.shape)
        self.method_quality_weight = np.where(
            self.method_valid, method_quality_weight, np.nan
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
        if self.count_source == "unsubtracted":
            del (
                self.wic_unsubtracted_counts,
                self.si12_unsubtracted_counts,
                self.si13_unsubtracted_counts,
                self.wic_unsubtracted_valid,
                self.si12_unsubtracted_valid,
                self.si13_unsubtracted_valid,
            )

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
            nc.count_source = self.count_source
            nc.method_quality_weight_method = self.method_quality_weight_method
            nc.proton_flux_source = self.proton_flux_source
            nc.proton_energy_model = self.proton_energy_model
            nc.proton_operation_order = self.proton_operation_order
            nc.proton_energy_uncertainty_method = (
                self.proton_energy_uncertainty_method
            )
            nc.proton_energy_coordinate_note = self.proton_energy_coordinate_note
            nc.count_uncertainty_mode = self.count_uncertainty_mode
            nc.count_uncertainty_method = self.count_uncertainty_method
            nc.spatial_smoothing_kernel = self.spatial_smoothing_kernel
            nc.wic_smoothing_width_pixels = self.wic_smoothing_width_pixels
            nc.si13_smoothing_width_pixels = self.si13_smoothing_width_pixels
            nc.spatial_smoothing_width_definition = (
                self.spatial_smoothing_width_definition
            )
            nc.spatial_smoothing_operation_order = (
                self.spatial_smoothing_operation_order
            )
            nc.spatial_smoothing_variance_method = (
                self.spatial_smoothing_variance_method
            )
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
                "si12": (self.si12, "counts"),
                "dsi12": (self.dsi12, "counts"),
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
            if self.spatial_smoothing_kernel != "none":
                fields.update({
                    "wic_smoothed": (self.wic_smoothed, "counts"),
                    "dwic_smoothed": (self.dwic_smoothed, "counts"),
                    "si13_smoothed": (self.si13_smoothed, "counts"),
                    "dsi13_smoothed": (self.dsi13_smoothed, "counts"),
                })
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
