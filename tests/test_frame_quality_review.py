from pathlib import Path
import sys

import numpy as np
import pandas as pd
import xarray as xr


SCRIPT_DIRECTORY = (
    Path(__file__).resolve().parents[1]
    / "scripts/background_quality_validation"
)
sys.path.insert(0, str(SCRIPT_DIRECTORY))

from annotate_frame_quality_review import (
    annotation_key,
    find_clicked_panel,
    load_annotations,
    load_completed,
    save_annotations,
    save_completed,
    set_manual_decision
)
from generate_frame_quality_review import (
    MANIFEST_COLUMNS,
    make_contact_sheet,
    rebuild_combined_manifest,
    select_edge_frames
)


def synthetic_frames(frame_count=5):
    shape = (frame_count, 4, 4)
    dates = np.datetime64("2000-05-18T00:00:00") + np.arange(frame_count) * np.timedelta64(2, "m")
    image = np.arange(np.prod(shape), dtype=float).reshape(shape) + 1

    return xr.Dataset(
        data_vars={
            "img": (("date", "row", "col"), image),
            "mlat": (("date", "row", "col"), np.full(shape, 70.0)),
            "sza": (("date", "row", "col"), np.full(shape, 100.0)),
            "dza": (("date", "row", "col"), np.full(shape, 30.0)),
            "frame_quality": ("date", np.arange(frame_count) % 3),
            "median_row_relative_spread": ("date", np.linspace(0.01, 0.10, frame_count)),
            "turn_on_reference_correlation": ("date", np.linspace(0.9, 0.1, frame_count)),
            "hv_mcp": ("date", np.full(frame_count, 1200.0)),
            "hv_phos": ("date", np.full(frame_count, 4500.0))
        },
        coords={
            "date": dates,
            "row": np.arange(shape[1]),
            "col": np.arange(shape[2]),
            "source_index": ("date", np.arange(frame_count)),
            "source_file": ("date", [f"frame_{index:03d}.sav" for index in range(frame_count)]),
            "northern_frame_index": ("date", np.arange(frame_count))
        }
    )


def manifest_row(automatic_quality=2):
    return pd.Series({
        "sheet": "wic/or_0085.png",
        "orbit": 85,
        "sensor": "WIC",
        "source_file": "frame_000.sav",
        "source_index": 0,
        "northern_frame_index": 0,
        "time": "2000-05-18T00:00:00",
        "automatic_quality": automatic_quality
    })


def test_select_edge_frames_returns_unique_chronological_indices():
    assert np.array_equal(
        select_edge_frames(50),
        np.r_[np.arange(18), np.arange(32, 50)]
    )
    assert np.array_equal(select_edge_frames(20), np.arange(20))
    assert select_edge_frames(0).size == 0


def test_select_edge_frames_rejects_invalid_arguments():
    with np.testing.assert_raises(ValueError):
        select_edge_frames(-1)
    with np.testing.assert_raises(ValueError):
        select_edge_frames(10, panels_per_end=0)


def test_contact_sheet_and_manifest_preserve_frame_identity(tmp_path):
    output_file = tmp_path / "wic" / "or_0085.png"
    manifest = make_contact_sheet(
        synthetic_frames(), orbit=85, sensor="wic",
        output_file=output_file, dpi=40
    )

    assert output_file.exists()
    assert list(manifest.columns) == MANIFEST_COLUMNS
    assert len(manifest) == 5
    assert manifest.sheet.unique().tolist() == ["wic/or_0085.png"]
    assert manifest.source_file.tolist() == [f"frame_{index:03d}.sav" for index in range(5)]
    assert manifest.automatic_quality.tolist() == [0, 1, 2, 0, 1]
    assert manifest.time.is_unique
    assert np.all((manifest.panel_left >= 0) & (manifest.panel_right <= 1))
    assert np.all((manifest.panel_top >= 0) & (manifest.panel_bottom <= 1))


def test_rebuild_combined_manifest_collects_completed_sheets(tmp_path):
    for sensor, orbit in (("wic", 85), ("si13", 86)):
        output_file = tmp_path / sensor / f"or_{orbit:04d}.png"
        manifest = make_contact_sheet(
            synthetic_frames(2), orbit=orbit, sensor=sensor,
            output_file=output_file, dpi=30
        )
        manifest.to_csv(output_file.with_suffix(".csv"), index=False)

    combined = rebuild_combined_manifest(tmp_path)

    assert len(combined) == 4
    assert (tmp_path / "manifest.csv").exists()
    assert set(combined.sensor) == {"WIC", "SI13"}


def test_click_mapping_returns_the_containing_panel():
    sheet_rows = pd.DataFrame([
        {"panel_left": 0.1, "panel_top": 0.2, "panel_right": 0.4, "panel_bottom": 0.5},
        {"panel_left": 0.5, "panel_top": 0.2, "panel_right": 0.8, "panel_bottom": 0.5}
    ])

    assert find_clicked_panel(sheet_rows, 0.2, 0.3).name == 0
    assert find_clicked_panel(sheet_rows, 0.6, 0.3).name == 1
    assert find_clicked_panel(sheet_rows, 0.45, 0.3) is None


def test_only_manual_exceptions_are_retained():
    annotations = {}
    row = manifest_row(automatic_quality=2)

    assert set_manual_decision(annotations, row, manual_quality=0) == "override"
    assert annotations[annotation_key(row)]["manual_quality"] == 0

    assert set_manual_decision(annotations, row, manual_quality=2) == "cleared"
    assert annotations == {}

    assert set_manual_decision(annotations, row, uncertain=True) == "uncertain"
    assert annotations[annotation_key(row)]["manual_quality"] == ""


def test_annotation_and_completion_state_round_trip(tmp_path):
    row = manifest_row(automatic_quality=2)
    annotations = {}
    set_manual_decision(annotations, row, manual_quality=0)
    annotation_file = tmp_path / "annotations.csv"
    save_annotations(annotation_file, annotations)

    completed = {
        "wic/or_0085.png": {
            "sheet": "wic/or_0085.png", "orbit": 85,
            "sensor": "WIC", "reviewed_utc": "2026-09-25T12:00:00+00:00"
        }
    }
    completed_file = tmp_path / "completed.csv"
    save_completed(completed_file, completed)

    loaded_annotations = load_annotations(annotation_file)
    loaded_completed = load_completed(completed_file)
    assert loaded_annotations[annotation_key(row)]["manual_quality"] == "0"
    assert "wic/or_0085.png" in loaded_completed
