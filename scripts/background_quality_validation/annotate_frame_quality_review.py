"""Annotate frame-quality contact sheets without opening the raw FUV data.

Click a panel, then press 0, 1, or 2 to record a manual quality override, or
``u`` to mark it uncertain. Press ``d`` to clear the selected annotation,
``n`` to mark the complete sheet reviewed and advance, ``p`` to go back, and
``q`` to save and quit. Decisions and completed-sheet state are written after
every action.
"""

#%% Imports

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


#%% Saved schemas

REQUIRED_MANIFEST_COLUMNS = {
    "sheet", "panel_index", "panel_left", "panel_top",
    "panel_right", "panel_bottom", "orbit", "sensor",
    "source_file", "source_index", "northern_frame_index",
    "time", "automatic_quality"
}
ANNOTATION_COLUMNS = [
    "sheet", "orbit", "sensor", "source_file", "source_index",
    "northern_frame_index", "time", "automatic_quality",
    "manual_quality", "status", "annotated_utc"
]
COMPLETED_COLUMNS = ["sheet", "orbit", "sensor", "reviewed_utc"]
SENSOR_ORDER = {"WIC": 0, "SI12": 1, "SI13": 2}


#%% Durable CSV storage

def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_csv_rows(filename):
    """Read a small state CSV or return no rows when it does not exist."""

    filename = Path(filename)
    if not filename.exists():
        return []
    with filename.open(newline="") as file:
        return list(csv.DictReader(file))


def atomic_write_rows(filename, columns, rows):
    """Atomically replace one annotation-state CSV."""

    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    temporary = filename.with_name(filename.name + ".partial")
    with temporary.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(filename)


def annotation_key(row):
    """Return a stable identity independent of panel position."""

    return str(row["sensor"]), str(row["source_file"]), str(row["time"])


def load_annotations(filename):
    return {
        annotation_key(row): row
        for row in read_csv_rows(filename)
    }


def load_completed(filename):
    return {
        str(row["sheet"]): row
        for row in read_csv_rows(filename)
    }


def save_annotations(filename, annotations):
    rows = sorted(
        annotations.values(),
        key=lambda row: (
            SENSOR_ORDER.get(row["sensor"], 99),
            int(row["orbit"]), row["time"]
        )
    )
    atomic_write_rows(filename, ANNOTATION_COLUMNS, rows)


def save_completed(filename, completed):
    rows = sorted(
        completed.values(),
        key=lambda row: (
            SENSOR_ORDER.get(row["sensor"], 99), int(row["orbit"])
        )
    )
    atomic_write_rows(filename, COMPLETED_COLUMNS, rows)


def set_manual_decision(annotations, manifest_row, manual_quality=None,
                        uncertain=False):
    """Set, replace, or clear one exception against the automatic flag."""

    row = manifest_row.to_dict() if hasattr(manifest_row, "to_dict") else dict(manifest_row)
    key = annotation_key(row)
    automatic = int(row["automatic_quality"])

    if not uncertain and int(manual_quality) == automatic:
        annotations.pop(key, None)
        return "cleared"

    annotations[key] = {
        "sheet": str(row["sheet"]),
        "orbit": int(row["orbit"]),
        "sensor": str(row["sensor"]),
        "source_file": str(row["source_file"]),
        "source_index": int(row["source_index"]),
        "northern_frame_index": int(row["northern_frame_index"]),
        "time": str(row["time"]),
        "automatic_quality": automatic,
        "manual_quality": "" if uncertain else int(manual_quality),
        "status": "uncertain" if uncertain else "override",
        "annotated_utc": utc_now()
    }
    return annotations[key]["status"]


#%% Manifest and click mapping

def load_manifest(review_directory):
    """Load and validate the combined server-generated panel manifest."""

    filename = Path(review_directory) / "manifest.csv"
    manifest = pd.read_csv(filename)
    missing = REQUIRED_MANIFEST_COLUMNS.difference(manifest.columns)
    if missing:
        raise ValueError(f"{filename} is missing {sorted(missing)}")
    if manifest.empty:
        raise ValueError(f"{filename} contains no review panels")

    manifest["sensor_order"] = manifest.sensor.map(SENSOR_ORDER).fillna(99)
    manifest = manifest.sort_values(
        ["sensor_order", "orbit", "panel_index"]
    ).reset_index(drop=True)
    return manifest


def find_clicked_panel(sheet_rows, x_normalized, y_normalized):
    """Return the manifest row containing one normalized image position."""

    inside = sheet_rows[
        (sheet_rows.panel_left <= x_normalized)
        & (x_normalized <= sheet_rows.panel_right)
        & (sheet_rows.panel_top <= y_normalized)
        & (y_normalized <= sheet_rows.panel_bottom)
    ]
    if inside.empty:
        return None
    return inside.iloc[0]


def panel_rectangle(row, image_width, image_height, **kwargs):
    """Make a display rectangle from normalized manifest bounds."""

    left = float(row["panel_left"]) * image_width
    top = float(row["panel_top"]) * image_height
    width = (float(row["panel_right"]) - float(row["panel_left"])) * image_width
    height = (float(row["panel_bottom"]) - float(row["panel_top"])) * image_height
    return Rectangle((left, top), width, height, fill=False, **kwargs)


#%% Interactive review

class FrameQualityReviewer:
    """Review pre-rendered contact sheets and persist only exceptions."""

    def __init__(self, review_directory, annotations_file, completed_file,
                 revisit_completed=False):
        self.review_directory = Path(review_directory)
        self.annotations_file = Path(annotations_file)
        self.completed_file = Path(completed_file)
        self.manifest = load_manifest(self.review_directory)
        self.annotations = load_annotations(self.annotations_file)
        self.completed = load_completed(self.completed_file)

        sheets = self.manifest[["sheet", "sensor_order", "orbit"]].drop_duplicates()
        sheets = sheets.sort_values(["sensor_order", "orbit"])
        if not revisit_completed:
            sheets = sheets[~sheets.sheet.isin(self.completed)]
        self.sheets = sheets.sheet.astype(str).tolist()
        self.position = 0
        self.selected = None
        self.image = None
        self.sheet_rows = None

        self.figure, self.axis = plt.subplots(figsize=(12, 12))
        self.figure.canvas.mpl_connect("button_press_event", self.on_click)
        self.figure.canvas.mpl_connect("key_press_event", self.on_key)

    def current_sheet(self):
        return self.sheets[self.position]

    def draw(self, message=""):
        """Render the current sheet and its existing local decisions."""

        self.axis.clear()
        self.selected = None
        sheet = self.current_sheet()
        image_file = self.review_directory / sheet
        if not image_file.exists():
            raise FileNotFoundError(f"Missing contact sheet: {image_file}")
        self.image = mpimg.imread(image_file)
        self.sheet_rows = self.manifest[self.manifest.sheet == sheet].copy()
        self.axis.imshow(self.image)
        self.axis.set_axis_off()

        height, width = self.image.shape[:2]
        for _, row in self.sheet_rows.iterrows():
            annotation = self.annotations.get(annotation_key(row))
            if annotation is None:
                continue
            colour = "cyan" if annotation["status"] == "uncertain" else "magenta"
            rectangle = panel_rectangle(
                row, width, height, edgecolor=colour, linewidth=4
            )
            self.axis.add_patch(rectangle)
            label = "?" if annotation["status"] == "uncertain" else f"M{annotation['manual_quality']}"
            self.axis.text(
                rectangle.get_x() + 4, rectangle.get_y() + 14, label,
                color=colour, fontsize=10, fontweight="bold",
                bbox={"facecolor": "black", "alpha": 0.65, "pad": 1}
            )

        sensor = str(self.sheet_rows.sensor.iloc[0])
        orbit = int(self.sheet_rows.orbit.iloc[0])
        instruction = (
            "Click panel; 0/1/2 override; u uncertain; d clear; "
            "n reviewed+next; p previous; q quit"
        )
        title = (
            f"{sensor} orbit {orbit:04d} — sheet {self.position + 1}/"
            f"{len(self.sheets)}\n{instruction}"
        )
        if message:
            title += f"\n{message}"
        self.axis.set_title(title, fontsize=11)
        self.figure.canvas.draw_idle()

    def on_click(self, event):
        if event.inaxes is not self.axis or event.xdata is None or event.ydata is None:
            return
        height, width = self.image.shape[:2]
        x_normalized = (event.xdata + 0.5) / width
        y_normalized = (event.ydata + 0.5) / height
        selected = find_clicked_panel(
            self.sheet_rows, x_normalized, y_normalized
        )
        if selected is None:
            return

        self.draw()
        self.selected = selected
        rectangle = panel_rectangle(
            selected, width, height, edgecolor="yellow", linewidth=4
        )
        self.axis.add_patch(rectangle)
        self.figure.canvas.draw_idle()

    def save_manual_quality(self, quality):
        if self.selected is None:
            self.draw("Select a panel first")
            return
        result = set_manual_decision(
            self.annotations, self.selected, manual_quality=quality
        )
        save_annotations(self.annotations_file, self.annotations)
        self.draw(f"Manual quality {quality}: {result}")

    def save_uncertain(self):
        if self.selected is None:
            self.draw("Select a panel first")
            return
        set_manual_decision(self.annotations, self.selected, uncertain=True)
        save_annotations(self.annotations_file, self.annotations)
        self.draw("Marked uncertain")

    def clear_selected(self):
        if self.selected is None:
            self.draw("Select a panel first")
            return
        removed = self.annotations.pop(annotation_key(self.selected), None)
        save_annotations(self.annotations_file, self.annotations)
        self.draw("Annotation cleared" if removed else "No annotation to clear")

    def complete_and_advance(self):
        sheet = self.current_sheet()
        first = self.sheet_rows.iloc[0]
        self.completed[sheet] = {
            "sheet": sheet,
            "orbit": int(first.orbit),
            "sensor": str(first.sensor),
            "reviewed_utc": utc_now()
        }
        save_completed(self.completed_file, self.completed)
        if self.position + 1 >= len(self.sheets):
            plt.close(self.figure)
            print("All selected contact sheets are reviewed")
            return
        self.position += 1
        self.draw("Previous sheet marked reviewed")

    def on_key(self, event):
        if event.key in ("0", "1", "2"):
            self.save_manual_quality(int(event.key))
        elif event.key == "u":
            self.save_uncertain()
        elif event.key in ("d", "delete", "backspace"):
            self.clear_selected()
        elif event.key == "n":
            self.complete_and_advance()
        elif event.key == "p":
            if self.position > 0:
                self.position -= 1
                self.draw()
        elif event.key == "q":
            plt.close(self.figure)

    def run(self):
        if not self.sheets:
            print("No unreviewed contact sheets remain")
            plt.close(self.figure)
            return
        self.draw()
        plt.show()


#%% Command line

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--review-dir", type=Path, required=True,
        help="transferred directory containing manifest.csv and sensor PNGs"
    )
    parser.add_argument(
        "--annotations", type=Path,
        help="override CSV; defaults to REVIEW_DIR/frame_quality_annotations.csv"
    )
    parser.add_argument(
        "--completed", type=Path,
        help="completed-sheet CSV; defaults to REVIEW_DIR/frame_quality_completed.csv"
    )
    parser.add_argument(
        "--revisit-completed", action="store_true",
        help="include sheets already recorded as reviewed"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    review_directory = args.review_dir.expanduser().resolve()
    annotations = (
        args.annotations.expanduser().resolve()
        if args.annotations is not None
        else review_directory / "frame_quality_annotations.csv"
    )
    completed = (
        args.completed.expanduser().resolve()
        if args.completed is not None
        else review_directory / "frame_quality_completed.csv"
    )
    reviewer = FrameQualityReviewer(
        review_directory, annotations, completed,
        revisit_completed=args.revisit_completed
    )
    reviewer.run()


if __name__ == "__main__":
    main()
