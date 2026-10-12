"""Convert a precipitation product into Hall and Pedersen conductance."""

#%% Imports

from copy import deepcopy
from pathlib import Path

import numpy as np
from icphysics import robinson_conductance
from icreader import load as icload
from netCDF4 import Dataset, date2num

from .grids import grid_mlt
from .precipitationdetector import ELECTRON_METADATA_ATTRIBUTES


#%% Conductance product

class ConductanceImage:
    """Apply a conductance model to an existing precipitation product."""

    def __init__(
        self,
        precipitation,
        conductance_function=None,
        conductance_model="robinson",
        conductance_provenance=None):
        
        source_file = None
        if isinstance(precipitation, (str, Path)):
            source_file = str(precipitation)
            precipitation = icload(precipitation)

        if getattr(precipitation, "product_type", "precipitation") != "precipitation":
            raise ValueError("ConductanceImage requires a precipitation product")

        if conductance_function is None:
            conductance_function = robinson_conductance

        self.product_type = "conductance"
        self.schema_version = 2
        self.precipitation_method = precipitation.method
        if self.precipitation_method == "hardy":
            for name in ELECTRON_METADATA_ATTRIBUTES:
                setattr(self, name, getattr(precipitation, name))
            self.E0_response_range_flag = np.asarray(precipitation.E0_response_range_flag, dtype=np.int8).copy()
        self.proton_flux_source = precipitation.proton_flux_source
        self.proton_energy_model = precipitation.proton_energy_model
        self.proton_energy_uncertainty_method = (
            precipitation.proton_energy_uncertainty_method
        )
        self.proton_energy_coordinate_note = precipitation.proton_energy_coordinate_note
        self.proton_response_energy_min = precipitation.proton_response_energy_min
        self.proton_response_energy_max = precipitation.proton_response_energy_max
        if self.proton_energy_model == "constant":
            self.proton_energy_constant = precipitation.proton_energy_constant
            self.proton_energy_uncertainty_constant = (
                precipitation.proton_energy_uncertainty_constant
            )
        self.conductance_model = conductance_model
        self.source_precipitation = source_file

        self.time = np.asarray(precipitation.time, dtype=object).copy()
        self.ssalon = np.asarray(precipitation.ssalon, dtype=float).copy()
        self.grid = deepcopy(precipitation.grid)
        self.shape = precipitation.shape
        if self.ssalon.shape != (self.shape[0],):
            raise ValueError("ssalon must contain one value per frame")
        self.kp = np.asarray(precipitation.kp, dtype=float).copy()
        self.kp_interval_start = np.asarray(precipitation.kp_interval_start, dtype="datetime64[ms]").copy()
        self.kp_provenance = dict(precipitation.kp_provenance)
        self.precipitation_provenance = dict(precipitation.physics_provenance)

        # Product 3 carries the precipitation state and its existing weight.
        for name in (
            "Ep_model", "Ep", "dEp", "Fp", "dFp",
            "E0", "dE0", "Fe", "dFe", "varE0Fe", "w",
        ):
            values = np.asarray(getattr(precipitation, name), dtype=float)
            if values.shape != self.shape:
                raise ValueError(f"{name} must have shape {self.shape}")
            setattr(self, name, values.copy())

        self.Ep_clipping_flag = np.asarray(
            precipitation.Ep_clipping_flag, dtype=bool
        ).copy()
        if self.Ep_clipping_flag.shape != self.shape:
            raise ValueError(f"Ep_clipping_flag must have shape {self.shape}")

        result = conductance_function(self.E0, self.Fe, self.dE0, self.dFe, self.varE0Fe)
        if not isinstance(result, dict):
            raise TypeError("conductance_function must return a dictionary")

        for name in ("P", "H", "dP", "dH"):
            if name not in result:
                raise ValueError(f"conductance_function did not return {name}")
            values = np.asarray(result[name], dtype=float)
            if values.shape != self.shape:
                raise ValueError(f"{name} must have shape {self.shape}")
            setattr(self, name, values.copy())

        provenance = {"module": conductance_function.__module__,
                      "function": getattr(conductance_function, "__name__", type(conductance_function).__name__)}
        provenance.update(conductance_provenance or {})
        self.conductance_provenance = provenance

    #%% NetCDF output

    def to_nc(self, filename):
        """Save the precipitation state and resulting conductance."""

        with Dataset(filename, "w", format="NETCDF4") as nc:
            nc.createDimension("time", self.shape[0])
            nc.createDimension("dim1", self.shape[1])
            nc.createDimension("dim2", self.shape[2])

            nc.product_type = self.product_type
            nc.schema_version = self.schema_version
            nc.precipitation_method = self.precipitation_method
            if self.precipitation_method == "hardy":
                for name in ELECTRON_METADATA_ATTRIBUTES:
                    nc.setncattr(name, getattr(self, name))
            nc.proton_flux_source = self.proton_flux_source
            nc.proton_energy_model = self.proton_energy_model
            nc.proton_energy_uncertainty_method = self.proton_energy_uncertainty_method
            nc.proton_energy_coordinate_note = self.proton_energy_coordinate_note
            nc.proton_response_energy_min = self.proton_response_energy_min
            nc.proton_response_energy_max = self.proton_response_energy_max
            if self.proton_energy_model == "constant":
                nc.proton_energy_constant = self.proton_energy_constant
                nc.proton_energy_uncertainty_constant = (
                    self.proton_energy_uncertainty_constant
                )
            nc.conductance_model = self.conductance_model
            if self.source_precipitation is not None:
                nc.source_precipitation = self.source_precipitation

            for name, value in self.precipitation_provenance.items():
                nc.setncattr(f"precipitation_physics_{name}", value)
            for name, value in self.conductance_provenance.items():
                nc.setncattr(f"conductance_{name}", value)
            for name, value in self.kp_provenance.items():
                nc.setncattr(f"kp_{name}", value)

            time_units = "seconds since 2000-01-01 00:00:00"
            time = nc.createVariable("time", "f8", ("time",))
            time[:] = date2num(self.time.tolist(), time_units, calendar="standard")
            time.units = time_units
            time.calendar = "standard"
            time.time_zone = "UTC"

            kp = nc.createVariable("Kp", "f4", ("time",))
            kp[:] = self.kp
            kp.units = "1"

            kp_start = nc.createVariable("Kp_interval_start", "f8", ("time",))
            kp_start[:] = date2num(self.kp_interval_start.astype(object).tolist(),
                                   time_units, calendar="standard")
            kp_start.units = time_units
            kp_start.calendar = "standard"

            ssalon = nc.createVariable("ssalon", "f4", ("time",))
            ssalon[:] = self.ssalon
            ssalon.units = "degrees"

            fields = {"Ep_model": (self.Ep_model, "keV"),
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
                      "w": (self.w, "1")}

            for name, (data, units) in fields.items():
                variable = nc.createVariable(name, "f4", ("time", "dim1", "dim2"), zlib=True)
                variable[:] = data
                variable.units = units

            if self.precipitation_method == "hardy":
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

            coordinates = {"xi": (self.grid.xi, "radians"),
                           "eta": (self.grid.eta, "radians"),
                           "mlat": (self.grid.lat, "degrees"),
                           "mlt": (grid_mlt(self.grid), "hours")}
            for name, (data, units) in coordinates.items():
                variable = grid.createVariable(name, "f8", ("dim1", "dim2"), zlib=True)
                variable[:] = data
                variable.units = units
