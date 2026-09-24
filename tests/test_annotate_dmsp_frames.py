from pathlib import Path
import sys

import numpy as np
import pandas as pd
import xarray as xr


SCRIPT_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "scripts/ratio_relation_validation"
)
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import annotate_dmsp_frames as annotation


def write_crossing(path, sza, mlt=23.0, image_data=True):
    path.mkdir()
    count = len(sza)
    time = np.datetime64("2001-01-01T00:00:00", "ns")
    dmsp_time = time + np.arange(count).astype("timedelta64[s]")
    data = xr.Dataset({
        "image_frame": ("sample", np.zeros(count, dtype=np.int32)),
        "image_time": ("sample", np.full(count, time)),
        "dmsp_sat": ("sample", np.full(count, "f15")),
        "dmsp_time": ("sample", dmsp_time),
        "dmsp_mlat": ("sample", np.linspace(65, 70, count)),
        "dmsp_mlt": ("sample", np.full(count, mlt)),
        "detector_row": ("sample", np.arange(count, dtype=np.int16)),
        "detector_column": ("sample", np.arange(count, dtype=np.int16)),
        "img_wic": (
            "sample", np.full(count, 100.0 if image_data else np.nan)
        ),
        "img_si13": (
            "sample", np.full(count, 5.0 if image_data else np.nan)
        ),
        "wic_sza": ("sample", np.asarray(sza, dtype=float)),
    }, attrs={
        "product_type": "dmsp_image_crossings",
        "representation": "detector",
    })
    data.to_netcdf(path / "or_0001.nc")


class FakeCSImage:
    product_type = "precipitation_cs"
    representation = "cs"
    time = np.array(["2001-01-01T00:00:00"], dtype="datetime64[ns]")

    def __init__(self):
        self.closed = False

    def close(self):
        self.closed = True


def deferred_annotation():
    return {
        "orbit": 1,
        "frame_id": 0,
        "img_time": "2001-01-01T00:00:00",
        "accepted": 2,
        "satellites": "F15",
        "n_dmsp_samples": 2,
        "annotated_utc": "2001-01-02T00:00:00+00:00",
    }


def test_minimum_sza_skips_without_saving_or_opening_an_image(tmp_path):
    matches = tmp_path / "matches"
    write_crossing(matches, [100.0, 104.9])
    output = tmp_path / "annotations.csv"

    summary = annotation.annotate(
        matches, tmp_path / "missing_images", output, minimum_sza=105
    )

    assert summary["non_dark_skipped"] == 1
    assert summary["presented"] == 0
    assert not output.exists()


def test_dark_candidate_uses_cs_context_and_saves_detector_key(
    tmp_path, monkeypatch
):
    matches = tmp_path / "matches"
    write_crossing(matches, [104.9, 105.0])
    images = tmp_path / "images"
    images.mkdir()
    (images / "or_0001.nc").touch()
    output = tmp_path / "annotations.csv"
    image = FakeCSImage()
    plotted = []

    monkeypatch.setattr(annotation, "open_cs_image", lambda filename: image)
    monkeypatch.setattr(
        annotation,
        "plot_frame",
        lambda *args, **kwargs: plotted.append((args, kwargs)) or object(),
    )
    monkeypatch.setattr(annotation, "read_keypress", lambda figure: "1")

    summary = annotation.annotate(
        matches, images, output, minimum_sza=105
    )
    saved = pd.read_csv(output)

    assert summary["presented"] == 1
    assert len(plotted) == 1
    assert image.closed
    assert saved.loc[0, "orbit"] == 1
    assert saved.loc[0, "frame_id"] == 0
    assert saved.loc[0, "img_time"] == "2001-01-01T00:00:00"
    assert saved.loc[0, "satellites"] == "F15"
    assert saved.loc[0, "accepted"] == 1


def test_without_sza_filter_preserves_legacy_mlt_deferred_status(tmp_path):
    matches = tmp_path / "matches"
    write_crossing(matches, [80.0, 81.0], mlt=12.0)
    output = tmp_path / "annotations.csv"

    summary = annotation.annotate(matches, tmp_path / "images", output)
    saved = pd.read_csv(output)

    assert summary["mlt_deferred"] == 1
    assert summary["presented"] == 0
    assert saved.loc[0, "accepted"] == 2


def test_deferred_candidate_remains_completed_by_default(tmp_path):
    matches = tmp_path / "matches"
    write_crossing(matches, [105.0, 106.0])
    output = tmp_path / "annotations.csv"
    annotation.save_annotation(output, deferred_annotation())

    summary = annotation.annotate(
        matches, tmp_path / "missing_images", output, minimum_sza=105
    )
    saved = pd.read_csv(output)

    assert summary["completed"] == 1
    assert summary["presented"] == 0
    assert saved["accepted"].tolist() == [2]


def test_revisit_deferred_replaces_row_without_duplicate(tmp_path, monkeypatch):
    matches = tmp_path / "matches"
    write_crossing(matches, [105.0, 106.0])
    images = tmp_path / "images"
    images.mkdir()
    (images / "or_0001.nc").touch()
    output = tmp_path / "annotations.csv"
    annotation.save_annotation(output, deferred_annotation())
    image = FakeCSImage()

    monkeypatch.setattr(annotation, "open_cs_image", lambda filename: image)
    monkeypatch.setattr(annotation, "plot_frame", lambda *args, **kwargs: object())
    monkeypatch.setattr(annotation, "read_keypress", lambda figure: "1")

    summary = annotation.annotate(
        matches, images, output, minimum_sza=105, revisit_deferred=True
    )
    saved = pd.read_csv(output)

    assert summary["presented"] == 1
    assert summary["deferred_replaced"] == 1
    assert len(saved) == 1
    assert saved.loc[0, "accepted"] == 1
    assert saved.loc[0, "annotated_utc"] != "2001-01-02T00:00:00+00:00"


def test_revisit_deferred_keeps_non_dark_row_untouched(tmp_path):
    matches = tmp_path / "matches"
    write_crossing(matches, [100.0, 104.9])
    output = tmp_path / "annotations.csv"
    annotation.save_annotation(output, deferred_annotation())

    summary = annotation.annotate(
        matches, tmp_path / "missing_images", output,
        minimum_sza=105, revisit_deferred=True,
    )
    saved = pd.read_csv(output)

    assert summary["non_dark_skipped"] == 1
    assert summary["deferred_replaced"] == 0
    assert len(saved) == 1
    assert saved.loc[0, "accepted"] == 2
    assert saved.loc[0, "annotated_utc"] == "2001-01-02T00:00:00+00:00"


def test_revisit_deferred_keeps_unusable_row_untouched(tmp_path):
    matches = tmp_path / "matches"
    write_crossing(matches, [105.0, 106.0], image_data=False)
    output = tmp_path / "annotations.csv"
    annotation.save_annotation(output, deferred_annotation())

    summary = annotation.annotate(
        matches, tmp_path / "missing_images", output,
        minimum_sza=105, revisit_deferred=True,
    )
    saved = pd.read_csv(output)

    assert summary["unusable_deferred_skipped"] == 1
    assert summary["detector_pair_missing_skipped"] == 0
    assert summary["deferred_replaced"] == 0
    assert len(saved) == 1
    assert saved.loc[0, "accepted"] == 2
    assert saved.loc[0, "annotated_utc"] == "2001-01-02T00:00:00+00:00"


def test_unannotated_missing_detector_pair_is_reported_without_decision(tmp_path):
    matches = tmp_path / "matches"
    write_crossing(matches, [105.0, 106.0], image_data=False)
    output = tmp_path / "annotations.csv"

    summary = annotation.annotate(
        matches, tmp_path / "missing_images", output, minimum_sza=105
    )

    assert summary["detector_pair_missing_skipped"] == 1
    assert summary["unusable_deferred_skipped"] == 0
    assert not output.exists()


def test_revisit_deferred_keeps_rejected_and_accepted_keys_completed(tmp_path):
    output = tmp_path / "annotations.csv"
    for orbit, accepted in ((1, 0), (2, 1), (3, 2)):
        row = deferred_annotation()
        row["orbit"] = orbit
        row["accepted"] = accepted
        annotation.save_annotation(output, row)

    completed = annotation.read_completed(output, revisit_deferred=True)

    assert completed == {
        (1, "2001-01-01T00:00:00", "F15"),
        (2, "2001-01-01T00:00:00", "F15"),
    }


def test_cs_frame_uses_absolute_time_if_indices_differ():
    image = FakeCSImage()
    image.time = np.array([
        "2000-12-31T23:59:00",
        "2001-01-01T00:00:00",
    ], dtype="datetime64[ns]")

    frame = annotation.cs_frame_for_time(
        image, np.datetime64("2001-01-01T00:00:00"), suggested_frame=0
    )

    assert frame == 1
