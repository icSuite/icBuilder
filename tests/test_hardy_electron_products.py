"""Hardy electron retrieval support, conditional errors and product contracts."""

import numpy as np
import pytest
from icphysics import hardy_electron_precipitation
from icreader import load, open_product
from netCDF4 import Dataset

from icbuilder.precipitationdetector import ELECTRON_METADATA_ATTRIBUTES, PrecipitationDetector
from icbuilder.conductancedetector import ConductanceDetector
from icbuilder.detectorcs import build_detector_cs_products
from icbuilder.precipitationimage import PrecipitationImage
from icbuilder.conductanceimage import ConductanceImage
from test_precipitationdetector import write_fuv_detector, kp_series, ORBIT_SCRIPT
from test_precipitationimage import _Binned, _kp_series, BASE_TIME
from icbuilder.grids import make_image_grids


@pytest.mark.parametrize("proton_model", ["constant", "hardy"])
@pytest.mark.parametrize("count_source", ["background_subtracted", "unsubtracted"])
def test_hardy_uses_wic_si12_when_si13_missing(tmp_path, proton_model, count_source):
    source = tmp_path / "fuv.nc"
    write_fuv_detector(source, two_frames=True)
    with Dataset(source, "r+") as nc:
        nc["wic_unsubtracted_counts"][:] = 10000
    kwargs = dict(method="hardy", kp_series=kp_series(), proton_energy_model=proton_model,
                  count_source=count_source)
    complete = PrecipitationDetector(source, **kwargs)
    # Invalidate all SI13 observations, errors and weights, including unsubtracted support.
    with Dataset(source, "r+") as nc:
        for name in ("si13_valid", "si13_unsubtracted_valid"):
            nc[name][:] = 0
        for name in ("si13_counts", "si13_variance", "si13_quality_weight"):
            nc[name][:] = np.nan
    missing = PrecipitationDetector(source, **kwargs)
    assert missing.method_valid[1].any()  # unmatched second SI13 frame remains useful
    for name in ("E0", "dE0", "Fe", "dFe", "varE0Fe", "Fp", "method_valid", "method_quality_weight"):
        np.testing.assert_allclose(getattr(missing, name), getattr(complete, name), equal_nan=True)
    expected = hardy_electron_precipitation(missing.kp[:, None, None], missing.mlt, missing.mlat)
    np.testing.assert_allclose(missing.E0[missing.method_valid], expected["mean_energy"][missing.method_valid])
    np.testing.assert_allclose(missing.dE0[missing.method_valid], 0)
    np.testing.assert_allclose(missing.varE0Fe[missing.method_valid], 0)


def test_hardy_missing_measurement_error_preserves_nominal(tmp_path):
    source = tmp_path / "fuv.nc"
    write_fuv_detector(source)
    nominal = PrecipitationDetector(source, method="hardy", kp_series=kp_series(), proton_energy_model="constant")
    with Dataset(source, "r+") as nc:
        nc["wic_variance"][:] = np.nan
    missing = PrecipitationDetector(source, method="hardy", kp_series=kp_series(), proton_energy_model="constant")
    for name in ("E0", "Fe", "method_valid"):
        np.testing.assert_allclose(getattr(nominal, name), getattr(missing, name), equal_nan=True)
    assert np.isnan(missing.dFe).all()
    np.testing.assert_allclose(missing.dE0[missing.method_valid], 0)
    filename = tmp_path / "precip.nc"
    missing.to_nc(filename)
    conductance = ConductanceDetector(filename)
    assert conductance.conductance_valid.any()
    assert not conductance.conductance_uncertainty_valid.any()


def test_hardy_restart_detects_electron_selection_and_provenance(tmp_path):
    source = tmp_path / "fuv.nc"
    write_fuv_detector(source)
    product = PrecipitationDetector(source, method="hardy", kp_series=kp_series(), proton_energy_model="constant", smoothed=False)
    output = tmp_path / "precip.nc"
    ORBIT_SCRIPT.save_precipitation_detector(product, output)
    kwargs = dict(filename=output, source_fuv_detector=source, proton_energy_model="constant",
                  proton_energy=2, proton_energy_uncertainty=0, smoothed=False)
    assert ORBIT_SCRIPT.precipitation_detector_file_status(**kwargs, method="hardy") == "complete"
    assert ORBIT_SCRIPT.precipitation_detector_file_status(**kwargs, method="image_ratio") == "mismatch"
    with Dataset(output, "r+") as nc:
        nc.electron_coefficient_source = "different-table"
    assert ORBIT_SCRIPT.precipitation_detector_file_status(**kwargs, method="hardy") == "mismatch"
    with Dataset(output, "r+") as nc:
        nc.electron_coefficient_source = product.electron_coefficient_source
        nc.renameVariable("E0_response_range_flag", "wrong_shape_response_flag")
        nc.createVariable("E0_response_range_flag", "i1", ("time",))[:] = 0
    assert ORBIT_SCRIPT.precipitation_detector_file_status(**kwargs, method="hardy") == "invalid"


def test_hardy_detector_conductance_and_cs_roundtrip(tmp_path):
    source = tmp_path / "fuv.nc"
    write_fuv_detector(source)
    p2 = PrecipitationDetector(source, method="hardy", kp_series=kp_series(), proton_energy_model="constant")
    p2_path = tmp_path / "precip.nc"
    p2.to_nc(p2_path)
    # Add a flagged detector pixel before Product 3 captures source identity.
    with Dataset(p2_path, "r+") as nc:
        flags = nc["E0_response_range_flag"][:]
        flags[p2.method_valid] = 2
        nc["E0_response_range_flag"][:] = flags
    p3 = ConductanceDetector(p2_path)
    p3_path = tmp_path / "conduct.nc"
    p3.to_nc(p3_path)
    assert p3.conductance_uncertainty_valid.any()
    pcs, ccs = build_detector_cs_products(p2_path, p3_path)
    for product, name, valid_name in ((pcs, "pcs", "method_valid"), (ccs, "ccs", "conductance_valid")):
        path = tmp_path / f"{name}.nc"
        product.to_nc(path)
        valid = getattr(product, valid_name)
        assert valid.any()
        np.testing.assert_allclose(product.E0_response_range_fraction[valid], 1)
        with open_product(path) as loaded:
            assert getattr(loaded, "method", getattr(loaded, "precipitation_method", None)) == "hardy"
            for attr in ELECTRON_METADATA_ATTRIBUTES:
                assert loaded.attrs[attr] == getattr(p2, attr)
            np.testing.assert_array_equal(loaded.read("E0_response_range_any"), product.E0_response_range_any)
            np.testing.assert_allclose(loaded.read("E0_response_range_fraction"), product.E0_response_range_fraction, equal_nan=True)


def test_hardy_bin_first_optional_si13_and_conductance_roundtrip(tmp_path):
    grid, _ = make_image_grids()
    wic = _Binned("WIC", [BASE_TIME], grid, 1000)
    si12 = _Binned("SI12", [BASE_TIME], grid, 4)
    si13 = _Binned("SI13", [BASE_TIME], grid, 6)
    kwargs = dict(wic=wic, si12=si12, method="hardy", kp_series=_kp_series(), proton_energy_model="constant")
    product = PrecipitationImage(**kwargs, si13=si13)
    missing = PrecipitationImage(**kwargs)
    for name in ("E0", "dE0", "Fe", "dFe", "varE0Fe", "w"):
        np.testing.assert_allclose(getattr(product, name), getattr(missing, name), equal_nan=True)
    path = tmp_path / "binned_precip.nc"
    missing.to_nc(path)
    restored = load(path)
    assert restored.method == "hardy"
    np.testing.assert_array_equal(restored.E0_response_range_flag, missing.E0_response_range_flag)
    conductance = ConductanceImage(path)
    path3 = tmp_path / "binned_conduct.nc"
    conductance.to_nc(path3)
    restored3 = load(path3)
    assert restored3.precipitation_method == "hardy"
    assert restored3.electron_energy_uncertainty_method == missing.electron_energy_uncertainty_method
    np.testing.assert_array_equal(restored3.E0_response_range_flag, missing.E0_response_range_flag)


def _pipeline_script(filename):
    import importlib.util
    from pathlib import Path

    path = Path(__file__).parents[1] / "scripts" / "pipeline" / filename
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("stage", ["bin_conductance", "detector_conductance", "precipitation_cs", "conductance_cs"])
@pytest.mark.parametrize("corruption", ["coefficients", "missing_flag", "missing_metadata"])
def test_hardy_downstream_restart_rejects_corruption(tmp_path, stage, corruption):
    if stage == "bin_conductance":
        grid, _ = make_image_grids()
        wic = _Binned("WIC", [BASE_TIME], grid, 1000)
        si12 = _Binned("SI12", [BASE_TIME], grid, 4)
        p2 = PrecipitationImage(wic, si12, "hardy", kp_series=_kp_series(), proton_energy_model="constant")
        p2_path = tmp_path / "p2.nc"
        p2.to_nc(p2_path)
        output = tmp_path / "p3.nc"
        ConductanceImage(p2_path).to_nc(output)
        script = _pipeline_script("make_conductance_orbit_files.py")
        status = lambda: script.conductance_matches_precipitation(output, p2_path)
        expected_complete = True
        expected_invalid = False
        flag_name = "E0_response_range_flag"
    else:
        source = tmp_path / "fuv.nc"
        write_fuv_detector(source)
        p2_path = tmp_path / "p2.nc"
        PrecipitationDetector(source, method="hardy", kp_series=kp_series(), proton_energy_model="constant").to_nc(p2_path)
        p3_path = tmp_path / "p3.nc"
        ConductanceDetector(p2_path).to_nc(p3_path)
        if stage == "detector_conductance":
            output = p3_path
            script = _pipeline_script("make_conductance_detector_orbit_files.py")
            status = lambda: script.conductance_detector_file_status(output, p2_path)
            flag_name = "E0_response_range_flag"
        else:
            pcs, ccs = build_detector_cs_products(p2_path, p3_path)
            pcs_path = tmp_path / "pcs.nc"
            pcs.to_nc(pcs_path)
            ccs_path = tmp_path / "ccs.nc"
            ccs.to_nc(ccs_path)
            script = _pipeline_script("make_detector_cs_orbit_files.py")
            if stage == "precipitation_cs":
                output = pcs_path
                status = lambda: script.precipitation_cs_file_status(output, p2_path)
            else:
                output = ccs_path
                status = lambda: script.conductance_cs_file_status(output, p3_path)
            flag_name = "E0_response_range_fraction"
        expected_complete = "complete"
        expected_invalid = "mismatch" if corruption == "coefficients" else "invalid"
    assert status() == expected_complete
    with Dataset(output, "r+") as nc:
        if corruption == "coefficients":
            nc.electron_coefficient_source = "wrong-source"
        elif corruption == "missing_metadata":
            nc.delncattr("electron_kp_policy")
        else:
            nc.renameVariable(flag_name, "removed_response_flag")
    assert status() == expected_invalid
