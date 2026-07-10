# Changelog

## Unreleased (working tree, 2026-07-10 — status-reporting correctness + QC strictness fix)

Companion pass to the same-day logging/status-reporting fixes in KUL_FWT and KUL_NIS.
`KUL_VBG.sh`'s own `task_exec`/`task_exec_soft` and `prep_log` were already correct
(confirmed by re-reading the live code: `prep_log` is freshly timestamped per invocation
via `$d`, and both functions already use `PIPESTATUS` and `exit 1`/`return 1` properly) —
the real gaps were in specific spots that bypassed those functions, plus VBG's QC verdict
logic being overly strict on one metric.

### KUL_VBG.sh
- Parallel tissue-prior warp block (4 backgrounded `antsApplyTransforms` calls): `wait ${_pid}`
  discarded the return value entirely, so a crashed warp job was invisible and downstream
  `mrcalc` steps would silently process missing/corrupt output. Now checks each job's exit
  status and aborts with a clear error if any failed.
- QC step invocation: `2>&1 | tee -a ${prep_log} || true` masked a genuine QC-script crash
  (not just a QC verdict of FAIL). Replaced with a `PIPESTATUS`-based check that logs a clear
  WARNING if `KUL_VBG_QC.py` itself crashed, while still treating QC as non-fatal to the
  overall run (a QC crash shouldn't abort a run that otherwise produced valid output).
- Donor-sharpening/hybrid-donor fallback logic (Task 2.5/2.7/2.7b): the primary attempt's
  `fslmaths` command substitutions were run with `2>/dev/null`, silencing the actual error
  before the `||` fallback kicked in — so if *both* the primary attempt and the fallback
  failed, there was no diagnostic trail at all. Changed to `2>>${prep_log}` so the error is
  preserved without changing the fallback control flow.

### KUL_VBG_QC.py
Real test-run feedback: a well-filled, visually good inpainting result was reported as
overall FAIL because `dark_voxel_fraction` (0.1246) exceeded its 0.05 threshold, even though
every other metric passed.

- Consolidated 3 duplicated/inconsistent threshold dicts (`render_metrics_summary`'s local
  4-metric dict missing `stitch_score`, the module-level `_THRESHOLDS` 5-metric dict, and the
  plain-text report's local 5-metric dict) plus 2 sets of hardcoded magic-number literals
  (console prints, plain-text "Overall pass/fail" section — which only checked 3 of 5 metrics,
  silently never failing on `fill_grad_vs_native_ratio` or `stitch_score`) into one
  module-level `_ARTIFACT_THRESHOLDS` source of truth.
- Added a `warn_only` flag to that shared threshold definition, set for
  `dark_voxel_fraction` specifically: it's computed as the fraction of lesion-mask voxels
  below native tissue's P10 across *all* tissue types (not WM-only, despite the old
  docstring) — since P10 by definition puts ~10% of native tissue below it already, any
  lesion mask that includes CSF-like/sulcal/ventricular margins (common for large lesions)
  will show a nontrivial "dark" fraction even on a well-filled result.
- Both the HTML and plain-text report generators now split `all_fails` from a separate
  `warns` list: a `warn_only` metric exceeding its threshold is still displayed/flagged
  clearly (new `.artifact-warn`/orange styling, "informational — not scored as failure") but
  no longer flips the overall verdict to FAIL on its own. Overall verdict is now one of PASS,
  PASS (with warnings), or FAIL (new `.overall-warn` CSS class for the middle case).
- Verified by hand: `dark_voxel_fraction=0.1246` with all 4 other metrics passing now yields
  `PASS (with warnings)` instead of `FAIL`, matching the real test-run case that prompted
  this fix.

## Unreleased (working tree, 2026-07-09)

- `KUL_VBG_QC.py`: added a `hex6()` helper that expands shorthand 3/4-digit
  hex colors (e.g. `#444`) to their 6-digit form before passing them to
  matplotlib's `spine.set_edgecolor()`. Fixes a compatibility issue with
  older matplotlib versions that can't parse shorthand hex codes.
- `KUL_lesion_overlap.py`: marked executable (`chmod +x`), no content change.
