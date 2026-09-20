"""Public icBuilder interface with workflow modules loaded on demand."""

from importlib import import_module

__all__ = [
    "PreImage",
    "BinnedImage",
    "PrecipitationImage",
    "PrecipitationDetector",
    "ConductanceDetector",
    "ConductanceImage",
    "make_image_grids",
    "make_wic_grid",
    "load_gfz_kp",
    "match_gfz_kp",
    "SplineImage",
    "load_zhang_paxton_lookup",
    "confun"
]


_PUBLIC_OBJECTS = {
    "PreImage": (".preimage", "PreImage"),
    "BinnedImage": (".binnedimage", "BinnedImage"),
    "PrecipitationImage": (".precipitationimage", "PrecipitationImage"),
    "PrecipitationDetector": (".precipitationdetector", "PrecipitationDetector"),
    "ConductanceDetector": (".conductancedetector", "ConductanceDetector"),
    "ConductanceImage": (".conductanceimage", "ConductanceImage"),
    "make_image_grids": (".grids", "make_image_grids"),
    "make_wic_grid": (".grids", "make_wic_grid"),
    "load_gfz_kp": (".kp", "load_gfz_kp"),
    "match_gfz_kp": (".kp", "match_gfz_kp"),
    "SplineImage": (".splineimage", "SplineImage"),
    "load_zhang_paxton_lookup": (
        ".zhang_paxton_lookup", "load_zhang_paxton_lookup"
    ),
    "confun": (".imagesat_e0_eflux_estimates", "E0_eflux_propagated"),
}


def __getattr__(name):
    """Load a public workflow object only when it is requested."""

    if name not in _PUBLIC_OBJECTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, object_name = _PUBLIC_OBJECTS[name]
    value = getattr(import_module(module_name, __name__), object_name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
