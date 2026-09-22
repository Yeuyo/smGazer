# smGazer Documentation

Single-molecule tracking (SMT) analysis of Sox2 at a fluorescently labelled transcription locus, with per-frame transcriptional state assignment from two reporter channels.

---

## 1. What the script does

For each cell in a directory tree, the script:

1. **Groups TIFF files** into cells by parsing filenames, separating imaging stacks from bead calibration images and assigning each to a channel.
2. **Computes channel registration offsets** from multi-colour bead images, using the 488 nm (GFP) channel as the reference.
3. **Segments GFP and RFP spots frame by frame**, records their presence, their mean intensity, the GFP spot coordinates, and whether GFP and RFP colocalise. Colocalisation is used as the proxy for the transcriptionally **ON** state.
4. **Localises and tracks single molecules** in the 640 nm channel (referred to throughout as Sox2) using the MTT-derived routines.
5. **Measures proximity** of each Sox2 localisation to the (last known) GFP locus, per frame and split by ON/OFF state.
6. **Estimates residence time** by fitting the survival curve of track durations with one- and two-component exponentials, then applying a photobleaching correction.
7. **Writes** per-cell spreadsheets and MATLAB figures, plus one pooled dwell-time spreadsheet for the whole run.

The script is a **script, not a function** — it runs top to bottom on a hard-coded input directory and leaves its results in the workspace and on disk.

---

## 2. Requirements

**MATLAB toolboxes**

| Toolbox | Functions used |
|---|---|
| Image Processing | `imquantize`, `multithresh`, `bwlabel`, `imadjust`, `imshow`, `imtool` |
| Statistics and Machine Learning | `pdist`, `squareform` |
| Optimization | `lsqnonlin`, `optimset` |

**External code, expected on the path**

- `tiffread` — reads a TIFF stack into a struct array with a `.data` field per frame. Must be on the MATLAB path already; the script does not add it.
- `Batch_MTT_code/` — must contain `localizeParticles_ASH` and `buildTracks2_ASH`. Added via

  ```matlab
  addpath(genpath(['.' filesep 'Batch_MTT_code' filesep]));
  ```

  This is a **relative** path, so MATLAB's current working directory must be the folder containing `Batch_MTT_code` when the script is run, regardless of where the data lives.

**Writing `.xls`** — `writecell`/`writematrix` to the legacy `.xls` extension. This works on current MATLAB but the format is deprecated; consider `.xlsx`.

---

## 3. Input data layout and naming

The script searches `mainDir` recursively (`**/*.tif`) and infers everything from filenames and folder names. The conventions are strict.

### Channel assignment (by filename suffix)

| Filename ends with | Channel | Variable |
|---|---|---|
| `488 nm.tif` | GFP reporter | slot 1 |
| `561 nm.tif` | RFP reporter | slot 2 |
| `640.tif` | SMT / Sox2 | slot 3 |

Note the inconsistency: the first two include `nm`, the third does not.

### Beads vs. cell images

A filename containing `BEADS` (case-insensitive) is treated as a **bead calibration image**; anything else is a **cell imaging stack**. Both are needed — a cell with either set incomplete is skipped (line 157).

### Cell number

Taken from the `numb`-th underscore-delimited token of the filename; the first run of digits in that token is converted to a number:

```matlab
fileName = strsplit(allFiles(n).name, '_');
cellNumber = str2double(regexp(fileName{numb},'[0-9\.]+','match'){1});
```

`numb` is set manually at the top (currently `2`). **This must match your naming scheme**, and the same token must also be the name of the folder holding that cell's images — the folder is only recorded when its last path component equals `fileName{numb}` (line 119). If it never matches, `cellFolder` stays `"0"` and the cell is silently dropped.

### ON/OFF folder tag

If the *folder path* contains `ON` or `OFF`, that string is stored in `allCells(n).state`. Note this is a plain `contains` test on the whole path, so any parent folder containing those letters (e.g. `…\Unanalysed\…` does not, but `…\CONTROL\…` would) will match. **This is due to the signals between ON and OFF states are manually assigned.**

### Indexing caveat

Cells are stored as `allCells(cellNumber)`, so cell numbers must be positive integers. Gaps produce empty struct entries, which are skipped.

---

## 4. Running it

1. Put MATLAB's current folder at the location containing `Batch_MTT_code/`.
2. Edit `mainDir` and `numb` (line 5).
3. Set the detection and tolerance parameters (Section 5).
4. Run. Results appear in a `Data` subfolder inside each cell's folder, plus `Summary_dwell_time.xls` in `mainDir`.

The comment at line 33 — "Specific GFP/RFP to use after run the analysis at least once" — describes the intended workflow: run once with `plotGFP`/`plotRFP` on, inspect the numbered spots in the saved per-frame figures, then re-run with `GFP2Use`/`RFP2Use` populated to force the choice of spot. **See Section 9.2 before using this mode — it overwrites the normal GFP/RFP detection. This is designed for cases where the user intend to manual assign GFP/RFP position**

---

## 5. Parameter reference

### 5.1 Spot detection and colocalisation

| Parameter | Value | Meaning |
|---|---|---|
| `minGFP` / `maxGFP` | 4 / 50 | Accepted size range (pixels) of a connected component to count as a GFP spot |
| `minRFP` / `maxRFP` | 4 / 50 | Same, for RFP |
| `distTol` | 4 | Distance threshold **in pixels** for GFP–RFP colocalisation *and* for Sox2-at-locus |
| `minBeads` | 4 | Minimum pixels for a bead in the SMT bead image (applied to the SMT channel only) |
| `gfpMobTol` | 12 | Maximum displacement (pixels) of the GFP spot between detections; larger jumps are rejected as a different spot |

At `impars.PixelSize = 0.11 µm`, `distTol = 4` corresponds to **0.44 µm** — roughly the diffraction limit. See Section 10.

### 5.2 Intensity thresholding

`imquantize(img, multithresh(img, Rng)) >= Tgt` — the image is split into `Rng+1` intensity classes by multilevel Otsu, and everything at or above class `Tgt` is kept.

| Parameter | Value | Meaning |
|---|---|---|
| `gfpTgt` / `gfpRng` | 14 / 20 | GFP: 21 classes, keep classes ≥ 14 |
| `rfpTgt` / `rfpRng` | 16 / 20 | RFP: 21 classes, keep classes ≥ 16 |

Both are wrapped in `try/catch` with a fallback of `multithresh(img, 5) >= 5`. `multithresh` accepts at most 20 levels, so `Rng = 20` sits exactly on the limit. The fallback is a **substantially more permissive** threshold, and it is applied silently — a frame that falls back is not flagged anywhere.

A frame-level sanity check follows: if more than 1% of the image passes threshold, the mask is discarded (`gfpLoc = 0`), on the assumption that the frame is background-dominated.

### 5.3 Timing

| Parameter | Value | Meaning |
|---|---|---|
| `frameIntervals` | 1 (commented alt: 30) | Number of SMT frames acquired between successive GFP/RFP frames |
| `frameTime` | 0 (commented alt: 0.606) | Seconds lost cycling through the GFP/RFP channels |
| `ExposureTime` | 500 ms | SMT exposure; sets `impars.FrameRate` |

With the current values the interleaved-acquisition timestamp correction at line 464 is a no-op. **Set both to match your acquisition** if GFP/RFP frames are interleaved with SMT frames.

### 5.4 Plotting

| Parameter | Value | Effect |
|---|---|---|
| `plotGFP` / `plotRFP` | 0 / 1 | Save per-frame annotated images with numbered spots |
| `plotNumber` | 10 | Plot only frames with index `< plotNumber` |
| `plotStart` | 1 | First frame to plot |
| `plotGap` | 1 | Plot every *n*-th frame |
| `plotHeatMap` | 0 | 512×512 occupancy heat map of all localisations |
| `plotAllGFPs` | 0 | Overlay all recorded GFP coordinates on one frame |

### 5.5 Manual overrides

| Parameter | Meaning |
|---|---|
| `GFP2Use` / `RFP2Use` | Per-detection-event spot indices to force, consumed one per frame that has detections |
| `manualGFP` | `[x y]` seed position; on frame 1 only, GFP spots further than `gfpMobTol` from this are rejected |

### 5.6 Localisation, tracking and photobleaching

| Parameter | Value | Notes |
|---|---|---|
| `bleach_rate` | 6.37 | Used in the residence-time correction; units ambiguous (see 9.5) |
| `LocalizationError` | −6.5 | Passed as `locpars.errorRate`, i.e. 10^−6.5 |
| `EmissionWavelength` | 580 nm | Feeds `impars.wvlnth` and the PSF width |
| `NumDeflationLoops` | 0 | Non-zero implies labelling density is too high |
| `MaxExpectedD` | 0.2 µm²/s | Tracking search radius driver |
| `NumGapsAllowed` | 1 | Frames a trajectory may skip |

`impars` (pixel size 0.11 µm, NA 1.49, PSF scale 1.35), `locpars` (9-pixel detection box, 50 optimisation iterations, no SNR/precision/density filtering) and `trackpars` (`searchExpFac` 1.2, `statWin` 10, `maxComp` 3) are passed straight through to the MTT routines.

---

## 6. Pipeline walkthrough

### Stage 1 — File cataloguing (lines 87–143)

Builds `allCells`, a struct array indexed by cell number with fields `name` (3×1 string of cell stacks), `bead` (3×1 string of bead images), `cell`, `state`, `folder`.

### Stage 2 — Bead registration (lines 167–243)

- SMT beads: thresholded at `multithresh(·,5) >= 2`, components smaller than `minBeads` discarded, centroids kept.
- RFP and GFP beads: thresholded at `multithresh(·,2) >= 3`; **no size filter**.
- Offsets are per-bead differences from the GFP centroids, NaN rows removed, then averaged and rounded to give `smtCorrection` and `rfpCorrection` (integer pixel shifts, `[Δx Δy]`).

Beads are matched **by `bwlabel` index**, i.e. by column-major scan order, not by nearest-neighbour matching. See Section 9.1.

### Stage 3 — Per-frame reporter analysis (lines 244–442)

For each frame: threshold GFP and RFP, label components, keep those within the size range.

For each candidate GFP spot:
- reject it if it is more than `gfpMobTol` from the **last frame in which GFP was accepted** (or, on frame 1, from `manualGFP` if set);
- record presence, centroid and mean intensity;
- loop over RFP spots, apply `rfpCorrection`, and on the first RFP within `distTol` mark the frame as colocalised, record the RFP mean intensity, and stop.

If no GFP is accepted in a frame, the coordinates are carried forward from the last **colocalised** frame (lines 383–389), and a second pass (lines 555–565) fills remaining gaps from the immediately preceding frame. The effect is that columns 4–5 always hold a "last known locus position".

### Stage 4 — SMT localisation and tracking (lines 444–458)

The whole FISH/640 stack is cast to `double` and passed to `localizeParticles_ASH`, then `buildTracks2_ASH`. Localisations are shifted by `smtCorrection` into GFP coordinates. `trackedPar` holds per-trajectory positions in µm, frame indices and timestamps.

No background subtraction, flat-field correction or drift correction is applied before localisation.

### Stage 5 — Dwell time (lines 460–539)

Track duration = `TimeStamp(end) − TimeStamp(1) + ExposureTime/1000`. A 1 s-binned histogram is converted to a survival curve `cdf_n`, truncated where it drops below 1%, and fitted:

- **1-component:** `cdf ≈ x1 + x2·exp(−t/x3)` (includes a constant offset)
- **2-component:** `cdf ≈ r1·exp(−t/r2) + (1−r1)·exp(−t/r3)` (no offset, amplitudes constrained to sum to 1)

Photobleaching correction:

```
TrueR  = (bleach_rate − 1) / (bleach_rate/x(3) − 1)
TrueR2 = (bleach_rate − 1) / (bleach_rate/max(r(2),r(3)) − 1)
```

### Stage 6 — Sox2-at-locus statistics (lines 567–789)

`onFrame` is built from colocalised frames, filling in runs separated by exactly `frameIntervals`. `offFrame` is **every other frame**. For each subset, the script counts localisations within `distTol` of the locus and the number of frames containing at least one.

A separate persistence measure, `sox2AppearLong`, counts a localisation as "specific binding" if a within-tolerance localisation in the previous recorded frame lay within **1 pixel** of it. This is computed independently of the MTT trajectories.

---

## 7. Data structure reference

### `cellData` — one row per GFP/RFP frame

| Col | Contents |
|---|---|
| 1 | GFP spot accepted (0/1) |
| 2 | RFP spot present (0/1) |
| 3 | GFP–RFP colocalised (0/1) → "ON" |
| 4 | GFP spot x (pixels, GFP frame of reference) |
| 5 | GFP spot y |
| 6 | Mean GFP intensity over the spot |
| 7 | Mean RFP intensity — **only written when colocalised**, otherwise 0 |

### `extraData` — 20 summary values, written to `_Summary.xls`

| # | Label as written | Computed as |
|---|---|---|
| 1 | Total number of Sox2 detected | `length(tracks)` — total **localisations**, not molecules |
| 2 | Number of Sox2 on GFP spot during ON state | Σ localisations within `distTol` over ON frames |
| 3 | % of time GFP spot is covered by Sox2 during ON | ON frames with ≥1 near localisation / ON frames |
| 4 | % of Sox2 on GFP spot during ON | ON frames with ≥1 near localisation / total localisations |
| 5 | Number of Sox2 on last detected GFP spot during OFF | as #2, over OFF frames |
| 6 | Number of Sox2 on GFP spot when only GFP detected | frames with col1=1, col2=0 |
| 7 | Number of Sox2 on last GFP spot when only RFP detected | frames with col1=0, col2=1 |
| 8 | …when neither GFP nor RFP detected | `#5 − #6 − #7` (derived, not measured) |
| 9 | Number of Sox2 detected during ON state | all localisations in ON frames |
| 10 | Number of frames in ON state | |
| 11 | Number of Sox2 detected during OFF state | all localisations in OFF frames |
| 12 | Number of frames in OFF state | |
| 13 | Minimum distance of Sox2 to GFP | `min(minDist)` over all frames — see 9.6 |
| 14 | % of time GFP spot covered by Sox2 during OFF | |
| 15 | Mean minimum distance, ON state | **µm** |
| 16 | Mean minimum distance, OFF state | **µm** |
| 17 | Visiting frequency, ON | `#2 / #10 × 100` |
| 18 | Visiting frequency, OFF | `#5 / #12 × 100` |
| 19 | Number of frames Sox2 on GFP spot, ON | |
| 20 | Number of frames Sox2 on GFP spot, OFF | |

### `tracks`

`[frame, x, y]` — all localisations, shifted into the GFP frame of reference. One row per localisation.

### `finalData`

`FileName`, `Cell_State_Counts` (= `cellData`), `SMT_Spots` (= `tracks`), `Extra` (= `extraData`). Built in the workspace but **never saved to disk**.

---

## 8. Outputs

Written to `<cellFolder>/Data/`, prefixed with the GFP filename minus extension:

| File | Contents |
|---|---|
| `_Data.xls` | Per-frame `cellData` columns 1–5 with headers |
| `_Summary.xls` | The 20 labelled `extraData` values |
| `_dwell_time.xls` | Long dwell time (`TrueR2`), fraction (`r(1)`), short dwell time (`r(3)`) |
| `_max_gfp_timing.xls` | Frame indices at maximum GFP intensity |
| `_sox_timing.xls` | Frame indices with ≥1 Sox2 near the locus |
| `_signal_intensity.fig` | GFP (×3) and RFP mean intensity bars |
| `_signal_intensity_with_sox.fig` | As above with Sox2 occurrence overlaid |
| `_sox2_history.fig` | Localisations per frame |
| `_time_lapse.fig` | GFP presence, Sox2 presence, colocalisation markers |
| `_sox_number_frames.fig` | Sox2 count and "specific binding" count per frame |
| `_sox_motif.fig` | Cumulative Sox2 occurrence |
| `_GFP_Frame_N.fig`, `_RFP_Frame_N.fig` | Per-frame annotated images (when plotting enabled) |

And in `mainDir`: **`Summary_dwell_time.xls`** — long dwell time, fraction, short dwell time, one row per cell index. Rows for skipped cells are blank, and there is **no cell-identifier column**, so rows are only interpretable by position.

The `×3` multiplier applied to GFP intensity before plotting (lines 837, 866) is a display scaling with no stated basis. It does not affect `_max_gfp_timing.xls`, but it does make the two channels in those figures non-comparable.

---

## 9. Known issues

### 9.1 Bead loops use the wrong count if different beads number exists in different channels

Beads are paired across channels **by label index**. `bwlabel` numbers components in column-major order, which is usually but not necessarily consistent when the thresholds, and therefore the detected object sets, differ per channel. A nearest-neighbour pairing would be more robust. The `minBeads` size filter is also applied to the SMT channel only, so noise blobs can enter the GFP and RFP centroid sets.

### 9.2 `break` exits the frame loop in the `GFP2Use` branch (line 430)

In the manual-selection branch there is no inner spot loop, so this `break` exits `for m = 1 : length(stackGFP)`. **Analysis stops at the first colocalised frame.** All downstream statistics are then computed on a truncated `cellData`. The automatic branch is fine — there the `break` correctly exits the `for p` RFP loop.

Also in that branch, `gfpUseInd` advances once per frame with detections; if `GFP2Use` is shorter than the number of such frames, the script errors with an index-out-of-bounds.

### 9.3 Cells with fewer than two colocalised frames error out

```matlab
for m = 1 : length(idx) - 1
  ...
end
onFrame = [onFrame; idx(m + 1)];
```

If `idx` is empty or has one element, the loop body never runs and `m` retains a **stale value from an earlier loop** (MATLAB does not clear it), so `idx(m+1)` throws. A cell where GFP and RFP never colocalise will halt the whole run rather than being skipped.

### 9.4 Dead code with a probable axis swap (lines 239–249)

`smtNew` and `rfpNew` are computed and never used. Within them, `smtCorrection(1)` (an x/column offset) is subtracted from row indices and `(2)` from column indices — the axes look transposed. Since the arrays are discarded this has no effect on results, but it should be removed or fixed rather than left as a trap.

Note also that at line 241 `smtCorrection` is still an N×2 matrix, so `smtCorrection(1)` takes the first element by linear indexing.

### 9.5 Residence-time fitting

- **The correction formula is undocumented.** `(b−1)/(b/τ−1)` is not the common form `1/τ_true = 1/τ_obs − 1/τ_bleach`. It may be correct for how `bleach_rate` is defined here, but the definition is not recorded anywhere — the comment says "unit in second" while the name says *rate*. The value `6.37` is hard-coded with no provenance, no measurement date and no per-dataset re-measurement.
- **The two fits use different model families.** The 1-component fit has a free constant offset `x(1)`; the 2-component fit does not and forces amplitudes to sum to 1. `TrueR` and `TrueR2` are therefore not directly comparable, and no model-selection test is performed between them.
- **Label/value mismatch when `r(3) > r(2)`.** `TrueR2` uses `max(r(2), r(3))` for the long component, but the spreadsheet always reports `r(3)` as "Short Dwell Time" and `r(1)` as "Fraction". If the optimiser returns `r(3) > r(2)`, the reported short time equals the long one and the fraction refers to the *other* component. Worth sorting the components explicitly and reporting the fraction that matches.
- **No fit diagnostics.** `lsqnonlin` is called with no bounds and its `exitflag`, residual and Jacobian are discarded, so negative or non-convergent parameters would pass through unnoticed. No confidence intervals are produced. Survival-curve bins are also not independent, so ordinary least squares on the CDF underestimates uncertainty; a maximum-likelihood fit to the durations would be better founded.
- **Fixed 1 s binning** with a 0.5 s exposure gives roughly two frames per bin, which coarsens exactly the short-lived population the two-component fit is meant to resolve.

### 9.6 `extraData(13)` is almost always zero

`minDist` is initialised to zeros and only written for frames containing localisations. `min(minDist)` over the whole vector therefore returns 0 whenever any frame lacks a detection. It should be `min(minDist(minDist > 0))` or restricted to frames present in `tracks(:,1)`.

### 9.7 `extraData(4)` mixes units

Labelled "Percentage of Sox2 on GFP spot during ON state", it divides a **frame count** (`orSox`) by a **localisation count** (`length(tracks)`). The result is not a percentage of anything interpretable. Compare with `extraData(2)`, which is a genuine localisation count.

Relatedly, `length(tracks)` on an N×3 matrix returns `max(N,3)`. With fewer than three localisations it silently returns 3.

### 9.8 `cellSummaryData` preallocated to the wrong size (line 810)

`cell(12,1)` is then filled to index 20. MATLAB grows it, so the output is correct, but the preallocation is misleading and hides the fact that 20 metrics exist.

---

## 10. Assumptions

**Locus position is carried forward.** During OFF frames the "GFP spot" is the last known position, so OFF-state proximity is measured against a stale coordinate. The longer the OFF stretch, the more chromatin motion accumulates and the more the measured distances are inflated. `gfpMobTol = 12 px` bounds the per-detection jump but not the drift across an undetected gap.

**Counts are localisations, not molecules.** `extraData(1)` and the ON/OFF Sox2 counts sum localisations across frames, so one molecule resident for 20 frames contributes 20. These are not counts of binding events and should not be compared across cells with different track-length distributions without normalising. The trajectory data needed to do this properly (`trackedPar`) is already computed but unused for these statistics.

**Two persistence measures, computed differently.** Dwell time comes from MTT trajectories; "specific binding" (`sox2AppearLong`) comes from a hand-rolled 1-pixel nearest-neighbour match between consecutive recorded frames, which bridges frames with no detections (the reset `soxInFrame = []` sits inside `if ~isempty(frame)`) and can increment more than once per frame. The two will not agree, and it is not documented which is authoritative. Deriving both from the trajectories would remove the discrepancy.