# Changelog

## Unreleased (working tree, 2026-08-14 — containerisation rebuilt for v2.0: Docker + Apptainer)

The v1.x `Docker/KUL_VBG_Dockerfile` could not produce a working VBG 2.0 image —
it installed FreeSurfer 6.0.0 (VBG 2.0 needs >= 7.3 for `mri_synthstrip`,
`mri_synthseg`, `segment_subregions` and `mri_segment_hypothalamic_subunits`),
still installed HD-BET (removed from VBG in Task 2.2), and built FastSurfer from
`env/fastsurfer.yml`, which is an empty placeholder upstream now. Rebuilt from
scratch, and extended to Apptainer for rootless multi-user servers.

### New — `Docker/`
- `Dockerfile` — multi-stage build on `ubuntu:24.04`: FreeSurfer 8.2.0, FSL
  6.0.7.23, ANTs 2.4.4.post20, MRtrix3 (CLI only), FastSurfer, KUL_VBG 2.0.
  Version pins mirror `KUL_Linux_setup/setup_environment.sh`, i.e. the versions
  the pipeline is validated against on bare metal. ANTs is built with gcc-12
  (its vendored ITKv5 snapshot does not compile under Ubuntu 24.04's default
  gcc-13 — ITK #4607), MRtrix3 with `-DMRTRIX_BUILD_GUI=OFF` (VBG calls no GUI
  tool), and FastSurfer's checkpoints are baked in at build time so nothing
  needs network or a writable image at run time.
- `KUL_VBG.def` + `build.sh --sif` — Apptainer SIF, produced *from* the Docker
  image so the two recipes cannot drift apart.
- `entrypoint.sh` — rewritten to never write inside the image (an Apptainer SIF
  is read-only; the old entrypoint `cp`-ed the licence into `$FREESURFER_HOME`
  and `/opt/FastSurfer`, which fails outright there). The FreeSurfer licence is
  now a bind mount at `/licence/license.txt`, reached through a symlink baked at
  `$FREESURFER_HOME/license.txt`. Also pins the tool environment explicitly,
  because Apptainer merges the host environment by default and a user's own
  `FSLDIR`/`FREESURFER_HOME` would otherwise shadow the container's.
- `fsl-packages-full.txt` — package manifest for `build.sh --fsl full`.
  Default is a measured subset (~1.5 GB vs ~11 GB): the complete set of FSL
  calls in this repo is `fslmaths`, `fslstats`, `fslreorient2std`, `fslswapdim`,
  `fslorient` and `imcp`, plus `MNI152_T1_1mm_brain.nii.gz` for QC. Installed as
  conda packages, so the shared-library closure is resolved rather than
  hand-copied.
- Removed `KUL_VBG_Dockerfile` and `Readme.txt` (superseded; in git history).

### Docker/Dockerfile — apply the official FreeSurfer 8.2.0 patches
Found by running `-M` in the container: every subregion segmentation
(hippocampal subfields, amygdala nuclei, thalamic nuclei, brainstem
substructures) silently failed, while Lausanne/Glasser/hypothalamic succeeded.

Cause is an upstream FreeSurfer packaging gap, not VBG: the 8.2.0 `.deb` ships
samseg's compiled gems bindings built only for Python 3.12
(`gemsbindings.cpython-312-*.so`), but FreeSurfer's bundled interpreter is
Python 3.8.13. `segment_subregions` therefore dies on import with
`ModuleNotFoundError: No module named 'samseg.gems.gemsbindings'`. Confirmed
with `dpkg -L freesurfer`, which lists only the cpython-312 build.

Upstream fixes this via `fs820_updates.sh`; the patch note for that file reads
"add gemsbindings.cpython python 3.8 library". The ubuntu24 patch set also
replaces `recon-all`, `mri_glmfit`, `rca-surfreg`, `segmentHA_T2.sh` and
others, so an unpatched install is an incomplete FreeSurfer, not a one-file
gap. The Dockerfile now applies it, pinned to `patch_data_20260625.tgz` for
reproducibility (`--build-arg FREESURFER_PATCH=` to bump).

Both the Dockerfile assertion and `build.sh`'s smoke/SIF checks now *run*
`segment_subregions --help` rather than only `command -v`-ing it — the binary
exists in a broken install, so the presence check passed and the image shipped
looking healthy. VBG itself behaved correctly throughout: it logged a non-fatal
WARNING and deliberately did not write `thalamic_nuclei.done`, so the work
retries on the next run.

### KUL_VBG_QC.py — crash in `render_metrics_summary` on a skipped stage
Found on the first full container run (sub-PT001, 2026-08-14). `render_metrics_summary`
filtered stages with `metrics["per_stage"][k]["p95_in_lesion"] > 0`, but
`compute_per_stage_metrics()` returns **every** field as `None` when a stage's
volume could not be loaded (its `vol is None` branch). Comparing `None > 0`
raises `TypeError` on Python 3, so the whole QC script died — after all six PNG
panels had been written, but *before* the text and HTML reports, leaving a QC
directory that looked almost complete but had no report in it.

On that run `Filled_MNI` was the skipped stage; the other 8 all produced values.

The console print, the HTML report and the plain-text report already guarded
this case (`if vol is not None` / `if m.get("p95_in_lesion") is None`);
`render_metrics_summary` was the only consumer that did not. Now guards for
`None` before the comparison, and defaults a `None` gradient to `0.0` for the
bar chart. Verified by re-running the function against the exact metrics JSON
that crashed, then re-running the whole script end to end: exit 0, verdict
`PASS (with warnings)`, both reports written.

Note this is not container-specific — it fires identically on a bare-metal run
whenever any stage volume is missing.

### KUL_VBG.sh — GPU detection fix in the `-P 3` hybrid block
Found while validating CPU-only container runs. The `-P 3` block ran
`nvidia-smi` unguarded, so on a machine without it `nvram` became an empty
string. The `[[ ! -z ${nvram} ]]` test correctly fell through to
`batch_fasu=4`, but the later `[ $nvram -lt 5500 ]` then evaluated as
`[ -lt 5500 ]` — which bash reports as "unary operator expected" and treats as
**false**, selecting the CUDA branch (`FaSu_cpu=""`) on a machine with no GPU.
Silent, and surfacing only as a confusing FastSurfer failure further on.

The `-P 2` block ~250 lines below already had the correct guard; `-P 3` had been
missed. Now uses the identical `command -v nvidia-smi || nvram=0` form, which
makes both the batch-size and the `--no_cuda` decisions correct. Verified by
hand: with `nvram=""` the old code selects CUDA, with `nvram=0` it selects
`--no_cuda`.

This is reachable outside containers too — any CPU-only machine running `-P 3` —
but running without `--nv`/`--gpus all` is a supported container mode, which is
what made it easy to hit.

### README.md
Dependency list corrected for v2.0 (FreeSurfer >= 7.3 required; HD-BET no longer
a dependency; Python packages listed), and the stale Docker section replaced
with Docker + Apptainer usage.

## Unreleased (working tree, 2026-07-11 — extend status-reporting pass to multiparc scripts)

Follow-up to the 2026-07-10 pass below, applying the same `.done`-marker audit to
`KUL_VBG_multiparc.sh` and the multiparc-related section of `KUL_VBG.sh` that were
flagged (but not checked) in that pass's final sweep.

Most of both files were already correct: `run`/`task_exec` (recon-all) is fatal
with `exit 1` on failure, so those `.done` touches were already unreachable after
a real failure, and the thalamic/brainstem/hippo-amygdala/hypothalamic subregion
steps already use `if run_soft/task_exec_soft; then touch ...; else WARNING; fi`
correctly. The one real gap in both files: the Lausanne2018 "no reference LUT for
this scale" fallback path does a raw `mv` (not routed through `run`/`task_exec`)
with nothing checking whether it succeeded before `multiscale_parc.done` gets
touched. Now checked explicitly in both `KUL_VBG_multiparc.sh` and `KUL_VBG.sh`.

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
