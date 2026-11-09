---
name: colony-count
description: |
  Count colonies in a colony formation assay (crystal violet) from photos of the dishes, with an
  interactive review: Claude proposes, the user removes/adds with a click, Claude recomputes.
  Outputs counts.csv: colony number + mean colony size + fold change vs CTRL. Figures only on request.
  TRIGGER: photos of dishes with violet colonies; "count the colonies"; "colony formation"; `/colony-count`.
  NOT for: crystal violet OD read in wells (plate reader), or Transwell membranes (use `/invasion-count`).
---

# colony-count

`SKILL_DIR` = the folder that contains this SKILL.md. Dependencies: `pip install numpy pandas pillow scipy scikit-image openpyxl` (R + tidyverse only for figures).

## Where it runs
- **Claude Code** (the user's own machine): everything as written below. `review` opens the page in their browser.
- **Cowork** (or any sandbox that cannot open the user's browser): work inside the folder the user shared, so they can see the files.
  - Install the dependencies first; if `pip install` refuses, add `--user` or `--break-system-packages`.
  - Do **not** run `review` (the page would be served inside the sandbox, out of the user's reach). Instead tell the user to open `<DIR>/review.html` from their folder with a double click.
  - On that page "Save" copies the corrections to the clipboard: the user pastes them in the chat, you write them to `<DIR>/corrections.json` and continue as after a save.
  - No R in the sandbox → skip figures unless `Rscript` exists.

## Steps
1. Create a results folder `<DIR>` (e.g. `colony_formation/<YYYY-MM-DD>/`) and copy the photos into it.
2. Look at each photo and write `<DIR>/design.csv` (`image,dish,line,condition,dose`). Dish order = rows from the top, then left to right. `dose` may be empty.
3. Ask the user only what the photo cannot tell: dose, cell line, which dish is the CTRL.
4. Run `detect` on every photo:
   ```
   python3 "$SKILL_DIR/colony_count.py" detect <photo> --design <DIR>/design.csv --out <DIR>
   ```
5. Read `<DIR>/<photo>_overlay.png` yourself first (green = colony, numbered per dish P1..Pn).
6. Show the overlay to the user and open the review page (run it in the background: it stays open until they press "Done"):
   ```
   python3 "$SKILL_DIR/colony_count.py" review <DIR>
   ```
7. Corrections, two ways:
   - Page: click a circle = remove/restore, click an empty spot = add. "Save and recompute" writes `<DIR>/corrections.json` and recomputes.
   - Chat: "P1 remove 12, 40; add one top left" → write `<DIR>/corrections.json` yourself, then run `apply`.
   - If the page was opened as a plain file (without `review`), "Save" downloads `corrections.json`: move it into `<DIR>/`.
8. Threshold wrong on the whole photo (too many / too few) → re-run `detect` with `--thr <n>` (default Otsu, min 25; higher = fewer colonies).
9. After saving (or after chat corrections) run `apply`, read `counts.csv` and always report **number, mean size and both fold changes**:
   ```
   python3 "$SKILL_DIR/colony_count.py" apply <DIR>
   ```
10. Figures only if the user asks: add `--plot --colour <treated colour>` (ggplot theme_bw, default `darkred`).
    Outputs `colonies_n.png`, `colonies_size.png`, `colonies_raw.png`, `colonies_fc.png`.
11. Repeat 5-9 until the user says ok.

## Files
- `corrections.json`: `{"<photo>": {"<dish>": {"remove": [id,...], "add": [[x,y,r],...]}}}` — pixel coordinates of the original photo.
- `counts.csv`: auto, removed, added, colonies, mean_size_px, fold_change, size_fold_change.
- `colonies.csv`: one row per colony with `size_px`.
- Size = colony pixels in the photo: comparable only within the same photo; the fold change is comparable across photos.
- Hand-added colonies: size = violet spot under the click (within 6 px), otherwise the click circle.
- Fold change = treated / CTRL of the same line (and same dose, if the column is filled in).

## Known limits
- Merged colonies in dense dishes: the watershed splits them only partly → review is mandatory.
- Pale, diffuse colonies: some stay under the threshold → try a lower `--thr`.
- The dish rim (outer 10% of the radius) is excluded: colonies on the rim must be added by hand.
- Auto-counting tends to miss very small colonies: always check the overlay before trusting the number.
- n=1 per condition: no statistics until there are biological replicates.
