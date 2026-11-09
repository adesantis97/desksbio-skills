---
name: invasion-count
description: |
  Measure an invasion/migration assay (Transwell, stained cells photographed under the microscope)
  from several fields per condition, with an interactive review: Claude proposes the threshold and
  one dot per cell, the user moves the threshold, clicks dots away or adds missing cells, drops bad fields.
  Primary measure = invaded cells per field, counted by shape (a large spread cell counts once).
  Second measure = % of the field covered by the stain (area). On confluent fields the count is a minimum.
  Outputs fields.csv (one row per photo), condition.csv (mean ± SD of the fields) and, over the
  whole tree, invasion_summary.csv + invasion.xlsx with fold change vs CTRL across biological replicates.
  TRIGGER: photos of stained Transwell membranes; "invasion assay"; "migration assay"; `/invasion-count`.
  NOT for: colonies on dishes (use `/colony-count`); crystal violet OD read in wells.
---

# invasion-count

`SKILL_DIR` = the folder that contains this SKILL.md. Dependencies: `pip install numpy pandas pillow scipy scikit-image openpyxl` (R + tidyverse only for figures).

## Where it runs
- **Claude Code** (the user's own machine): everything as written below. `review` opens the page in their browser.
- **Cowork** (or any sandbox that cannot open the user's browser): work inside the folder the user shared, so they can see the files.
  - Install the dependencies first; if `pip install` refuses, add `--user` or `--break-system-packages`.
  - Do **not** run `review` (the page would be served inside the sandbox, out of the user's reach). Instead tell the user to open `<ROOT>/<rep>/<line>/<cond>/_invasion/review.html` from their folder with a double click.
  - On that page "Save" copies the corrections to the clipboard: the user pastes them in the chat, you write them to `<ROOT>/<rep>/<line>/<cond>/_invasion/corrections.json` and continue as after a save.
  - No R in the sandbox → skip figures unless `Rscript` exists.

## The folder tree
```
<ROOT>/<replicate>/<line>/<condition>/*.bmp   ← the 20X photos, ~10 per condition
                                     _invasion/  ← written by the skill, inside the condition
```
Names are free: only the three levels matter. Photos at another magnification (`..._10X.bmp`) are
skipped automatically: use `--mag 10` if those are the ones you need. Files without a `NNX` tag are always used.

## Steps
1. Put the photo tree in the project folder.
2. Ask the user only what the names do not tell: which condition is the CTRL, which lines to compare.
3. Measure everything:
   ```
   python3 "$SKILL_DIR/invasion_count.py" run <ROOT>
   ```
4. **Open the review of two conditions first** (one dense, one sparse: calibrate on the pair):
   ```
   python3 "$SKILL_DIR/invasion_count.py" review <ROOT>/<rep>/<line>/<cond>
   ```
   The page: threshold slider (green = stained), one yellow dot per counted cell, cells and area at the top.
   - Move the slider, then "Save and recompute": the dots are recounted at the new threshold (the page reloads).
   - Click a dot = remove it (red ×); click a cell without a dot = add one (blue); "exclude" drops a field.
   - "Save and recompute" writes `_invasion/corrections.json` and recomputes.
   Dot clicks belong to one threshold: changing the threshold after clicking discards them.
5. Take the threshold the user approved and **re-run everything with it** (clicks stay only on the reviewed conditions):
   ```
   python3 "$SKILL_DIR/invasion_count.py" run <ROOT> --thr <threshold>
   ```
6. Aggregate and show the table (cells first, area beside):
   ```
   python3 "$SKILL_DIR/invasion_count.py" aggregate <ROOT> --ctrl CTRL
   ```
7. Figures only if the user asks: `--plot [--colour <colour>]`.
8. Look at 2-3 `_invasion/*_overlay.jpg` yourself before saying the number is good.

## The threshold is the only choice
- Within the same replicate **and** line the threshold must be **the same for every condition**:
  that is what makes counts and percentages comparable. Default 0.10, calibrated on 20X crystal-violet membranes.
- Across replicates it may change (different staining and exposure): the fold change vs the CTRL of
  the same replicate absorbs the difference. That is why aggregation is done on fold changes, not percentages.
- A wrong threshold shows in the overlay: green on the beige membrane = too low, holes inside cells = too high.

## Files
- `_invasion/fields.csv` — one row per photo: `cells`, `cells_auto`, `removed`, `added`, `area_pct`, `iod`, `used`, `note`, `threshold`.
- `_invasion/condition.csv` — `cells_mean`, `cells_sd`, `cells_cv_pct`, `confluent_fields`, then the same for the area.
- `_invasion/corrections.json` — `{"threshold": 0.12, "exclude": ["photo.bmp"], "cells_thr": 0.12, "cells": {"photo.bmp": {"remove": [3, 17], "add": [[x, y]]}}}`;
  `remove` = index in the dot list of `detect.json`, `add` = coordinates at width 1536. To correct from chat, write this file and run `apply`.
- `_invasion/<photo>_overlay.jpg` — mask in green and the dots in yellow over the photo. Always look at it.
- `<ROOT>/invasion_fields.csv · invasion_conditions.csv · invasion_summary.csv` — from `aggregate`, plus `invasion.xlsx` with the three as sheets (summary, conditions, fields): tell the user it is there.
- `area_pct` = % of the field above threshold, after removing objects < 20 px (debris, pores).
- `iod` = mean stain index over the field: independent of the threshold, a cross-check.
- `cells` = one dot per lobe of the stained mask (peaks of the distance transform, min distance = one cell radius,
  22 px at width 1536). Every photo is resized to width 1536 first: fields shot at 3072 and 1536 px are comparable.
- `fold_change` = mean cells of the condition / mean cells of the CTRL **of the same replicate and line**;
  `area_fold_change` = the same on the area.
- `p_vs_ctrl` = one-sample t-test of the fold changes against 1, across biological replicates (needs n ≥ 3).

## Known limits
- Confluent field (> 60% area): cells merge into a mesh, the count is a **minimum** (marked in `note` and
  `confluent_fields`). Say it next to the number; the area saturates too.
- Out-of-focus or very dark cells turn brownish and the G−B index catches them poorly: check the overlay.
- Air bubbles, folds and tears in the membrane: exclude them by hand on the page, do not correct them.
- Different number of fields per condition (9 vs 10): the mean holds, but write it in the report.
- Same field of view assumed for every photo (same objective and camera): if the magnification differs, do not compare counts.
- Write "invaded cells per field (shape count)", and the area as "invaded area (%)".
