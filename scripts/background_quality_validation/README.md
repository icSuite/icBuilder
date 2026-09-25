# FUV frame-quality review

This workflow checks whether the automatic `fuvpy` quality classifier accepted
bad frames near the beginning or end of an orbit. It does not modify Product 1
or the classifier.

The server step reads the raw IDL data serially and recalculates the same frame
quality used by the production background pipeline. It writes one 6-by-6
contact sheet per sensor/orbit containing the unique first and last 18 northern
frames. Borders show the automatic quality: red `0`, yellow `1`, green `2`.

```bash
mkdir -p ~/IMAGE_FUV/frame_quality_review
nohup python -u scripts/background_quality_validation/generate_frame_quality_review.py \
    --base-input ~/IMAGE_FUV \
    --output ~/IMAGE_FUV/frame_quality_review \
    > ~/IMAGE_FUV/frame_quality_review/generate.log 2>&1 < /dev/null &
```

The run is restartable. Existing orbit sheets with their per-sheet manifests
are skipped unless `--overwrite` is supplied. `failures.csv` records failed
sensor/orbit reads, while `manifest.csv` combines all completed sheets.

Transfer the review directory, except `generate.log` if it is large, to the
local machine. Annotation uses only the PNG files and CSV metadata; it never
opens the raw IDL or Product 1 files.

```bash
python scripts/background_quality_validation/annotate_frame_quality_review.py \
    --review-dir /path/to/frame_quality_review
```

Click a panel and use:

- `0`, `1`, or `2`: assign a manual quality;
- `u`: mark the frame uncertain;
- `d`: delete the selected decision;
- `n`: record the whole sheet as reviewed and advance;
- `p`: return to the previous sheet in the current session;
- `q`: save and quit.

`frame_quality_annotations.csv` stores only disagreements with the automatic
quality plus uncertain frames. A decision equal to the automatic quality is
removed because it is not an exception. `frame_quality_completed.csv` records
every completely reviewed sheet, including sheets with no exceptions. Both
files are updated after each action, so the review can be resumed safely.

Annotations use `(sensor, source_file, time)` as the stable frame identity.
They are validation evidence, not automatic production exclusions. After the
review, the exception table can be used to diagnose and revise the classifier
or, following an explicit scientific decision, to build an exclusion list.
