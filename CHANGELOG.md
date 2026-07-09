# Changelog

## Unreleased (working tree, 2026-07-09)

- `KUL_VBG_QC.py`: added a `hex6()` helper that expands shorthand 3/4-digit
  hex colors (e.g. `#444`) to their 6-digit form before passing them to
  matplotlib's `spine.set_edgecolor()`. Fixes a compatibility issue with
  older matplotlib versions that can't parse shorthand hex codes.
- `KUL_lesion_overlap.py`: marked executable (`chmod +x`), no content change.
