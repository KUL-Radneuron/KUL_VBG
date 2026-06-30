#!/usr/bin/env python3
"""
KUL_lesion_overlap.py

Compute per-label overlap between a parcellation image and a binary lesion mask.
Standalone — can be called directly or from KUL_VBG.sh via KUL_lesion_overlap_report().

Usage:
    python3 KUL_lesion_overlap.py \\
        --parc   parc.nii.gz \\
        --lesion lesion_mask.nii.gz \\
        --out    report.txt \\
        [--lut   lut.txt] \\
        [--name  "atlas label"] \\
        [--html  report.html]

LUT format (either accepted):
    2-column TSV:  label_id<TAB>label_name
    6-column FS:   label_id  label_name  R  G  B  A   (# comment lines skipped)

Output: tab-separated table, only rows where overlap > 0,
        sorted by %_of_lesion descending.
        When --html is given the section is appended to (or creates) an HTML report.

Exit codes:
    0  success
    1  input file missing
    2  image geometry mismatch (different voxel sizes or array shape)
"""

import argparse
import sys
import os
import re
import datetime
import numpy as np
import nibabel as nib
from nibabel.processing import resample_from_to


def load_lut(path):
    """Return dict {label_id: label_name}. Handles 2-col TSV and 6-col FS LUT."""
    lut = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                lut[int(parts[0])] = parts[1]
            except ValueError:
                continue
    return lut


def vox_vol_mm3(img):
    zooms = img.header.get_zooms()[:3]
    return float(zooms[0]) * float(zooms[1]) * float(zooms[2])


# ── HTML helpers ─────────────────────────────────────────────────────────────

_CSS = """
body{font-family:Arial,sans-serif;font-size:13px;margin:20px;background:#f8f9fa;color:#222}
h1{color:#2c3e50;border-bottom:2px solid #3498db;padding-bottom:6px}
.meta{background:#fff;border:1px solid #dee2e6;border-radius:4px;padding:10px 14px;
      margin-bottom:16px;font-size:12px;color:#555}
details{background:#fff;border:1px solid #dee2e6;border-radius:4px;margin-bottom:10px}
summary{cursor:pointer;padding:10px 14px;font-weight:bold;font-size:13px;
        background:#f1f3f5;border-radius:4px;list-style:none;user-select:none}
summary::before{content:"▶ ";font-size:10px;margin-right:4px}
details[open] summary::before{content:"▼ "}
.tbl-wrap{padding:10px 14px 14px}
table{border-collapse:collapse;width:100%;font-size:12px}
th{background:#3498db;color:#fff;padding:5px 8px;text-align:left;cursor:pointer;
   white-space:nowrap;user-select:none}
th:hover{background:#2980b9}
td{padding:4px 8px;border-bottom:1px solid #eee;white-space:nowrap}
tr:hover td{background:#f0f7ff}
.bar-cell{width:120px}
.bar{height:10px;background:#3498db;border-radius:2px;min-width:2px}
.badge{display:inline-block;background:#27ae60;color:#fff;border-radius:3px;
       padding:1px 6px;font-size:11px;margin-left:6px}
.lesion-info{font-size:11px;font-weight:normal;color:#555;margin-left:10px}
"""

_JS = """
function sortTable(tbl, col, asc) {
  const tbody = tbl.tBodies[0];
  const rows  = Array.from(tbody.rows);
  rows.sort((a, b) => {
    const va = a.cells[col].dataset.v ?? a.cells[col].textContent.trim();
    const vb = b.cells[col].dataset.v ?? b.cells[col].textContent.trim();
    const na = parseFloat(va), nb = parseFloat(vb);
    if (!isNaN(na) && !isNaN(nb)) return asc ? na - nb : nb - na;
    return asc ? va.localeCompare(vb) : vb.localeCompare(va);
  });
  rows.forEach(r => tbody.appendChild(r));
}
document.addEventListener('DOMContentLoaded', () => {
  document.querySelectorAll('table').forEach(tbl => {
    let sortCol = -1, sortAsc = false;
    tbl.querySelectorAll('th').forEach((th, i) => {
      th.addEventListener('click', () => {
        sortAsc = (sortCol === i) ? !sortAsc : false;
        sortCol = i;
        sortTable(tbl, i, sortAsc);
      });
    });
  });
});
"""

_SENTINEL = "<!-- KUL_SECTIONS_END -->"


def _html_header(subject="", les_vox=0, les_mm3=0.0):
    date_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    subj_line = f"<b>Subject:</b> {subject} &nbsp;|&nbsp;" if subject else ""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>KUL Lesion Overlap Report</title>
<style>{_CSS}</style>
</head>
<body>
<h1>KUL Lesion Overlap Report</h1>
<div class="meta">
  {subj_line}
  <b>Date:</b> {date_str} &nbsp;|&nbsp;
  <b>Lesion volume:</b> {les_vox:,} vox &nbsp;({les_mm3:.1f} mm³)
</div>
{_SENTINEL}
<script>{_JS}</script>
</body>
</html>
"""


def _html_section(name, rows, les_vox, les_mm3, parc_path, lut_path, open_tag=False):
    open_attr = " open" if open_tag else ""
    n = len(rows)
    max_pct = max((r[6] for r in rows), default=1.0) or 1.0

    def bar(pct):
        w = max(2, int(100 * pct / max_pct))
        return f'<div class="bar" style="width:{w}px"></div>'

    hdr_cols = ["Structure", "Lbl vox", "Lbl mm³", "Ovl vox", "Ovl mm³", "% of lbl", "% of lesion", ""]
    thead = "".join(f"<th>{c}</th>" for c in hdr_cols)

    trows = []
    for name_r, lbl_vox, lbl_mm3, ovl_vox, ovl_mm3, pct_lbl, pct_les in rows:
        trows.append(
            f"<tr>"
            f"<td>{name_r}</td>"
            f'<td data-v="{lbl_vox}">{lbl_vox:,}</td>'
            f'<td data-v="{lbl_mm3:.1f}">{lbl_mm3:.1f}</td>'
            f'<td data-v="{ovl_vox}">{ovl_vox:,}</td>'
            f'<td data-v="{ovl_mm3:.1f}">{ovl_mm3:.1f}</td>'
            f'<td data-v="{pct_lbl:.1f}">{pct_lbl:.1f}%</td>'
            f'<td data-v="{pct_les:.1f}">{pct_les:.1f}%</td>'
            f'<td class="bar-cell">{bar(pct_les)}</td>'
            f"</tr>"
        )

    tbody = "\n".join(trows)
    lesion_info = f'<span class="lesion-info">(lesion: {les_vox:,} vox | {les_mm3:.1f} mm³)</span>'
    badge = f'<span class="badge">{n}</span>'

    return (
        f'<details{open_attr}>\n'
        f'<summary>{name}{badge}{lesion_info}</summary>\n'
        f'<div class="tbl-wrap">\n'
        f'<table><thead><tr>{thead}</tr></thead>\n'
        f'<tbody>\n{tbody}\n</tbody></table>\n'
        f'<p style="font-size:11px;color:#888;margin-top:6px">'
        f'Parc: {parc_path} &nbsp;|&nbsp; LUT: {lut_path or "none"}</p>\n'
        f'</div>\n'
        f'</details>\n'
    )


def append_html_section(html_path, section_html, les_vox, les_mm3, subject=""):
    if not os.path.exists(html_path):
        with open(html_path, 'w') as fh:
            fh.write(_html_header(subject=subject, les_vox=les_vox, les_mm3=les_mm3))

    with open(html_path) as fh:
        content = fh.read()

    if _SENTINEL not in content:
        # fallback: insert before </body>
        content = content.replace("</body>", f"{_SENTINEL}\n</body>")

    content = content.replace(_SENTINEL, section_html + _SENTINEL)

    with open(html_path, 'w') as fh:
        fh.write(content)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Per-label lesion overlap report")
    ap.add_argument("--parc",         required=True, help="Parcellation NIfTI (.nii/.nii.gz)")
    ap.add_argument("--lesion",       required=True, help="Binary lesion mask NIfTI")
    ap.add_argument("--out",          required=True, help="Output report text file")
    ap.add_argument("--lut",          default=None,  help="Primary label LUT file (optional)")
    ap.add_argument("--lut-fallback", default=None,  help="Fallback LUT for labels absent from --lut (e.g. FreeSurferColorLUT.txt)")
    ap.add_argument("--name",         default="",    help="Atlas label printed in report header")
    ap.add_argument("--html",         default=None,  help="HTML report file (created/appended per atlas)")
    ap.add_argument("--subject",      default="",    help="Subject ID (written to HTML header on first call)")
    args = ap.parse_args()

    # ── Input checks ────────────────────────────────────────────────────────────
    for path in (args.parc, args.lesion):
        if not os.path.isfile(path):
            print(f"ERROR: file not found: {path}", file=sys.stderr)
            sys.exit(1)

    parc_img   = nib.load(args.parc)
    lesion_img = nib.load(args.lesion)

    if parc_img.shape[:3] != lesion_img.shape[:3]:
        print(
            f"WARNING: shape mismatch — parc {parc_img.shape[:3]} vs lesion {lesion_img.shape[:3]}, resampling lesion to parc space",
            file=sys.stderr,
        )
        lesion_img = resample_from_to(lesion_img, parc_img, order=0)

    parc_vox   = np.round(parc_img.header.get_zooms()[:3], 4)
    lesion_vox = np.round(lesion_img.header.get_zooms()[:3], 4)
    if not np.allclose(parc_vox, lesion_vox, atol=0.01):
        print(
            f"WARNING: voxel size mismatch after resampling — parc {parc_vox} vs lesion {lesion_vox}, proceeding anyway",
            file=sys.stderr,
        )

    # ── Load data ────────────────────────────────────────────────────────────────
    parc   = np.asarray(parc_img.dataobj, dtype=np.int32)
    lesion = np.asarray(lesion_img.dataobj, dtype=np.float32) > 0.5
    vox_mm3 = vox_vol_mm3(parc_img)

    lut = load_lut(args.lut) if args.lut else {}
    if args.lut_fallback:
        fallback = load_lut(args.lut_fallback)
        # only fill in labels not already in the primary LUT
        for k, v in fallback.items():
            if k not in lut:
                lut[k] = v

    les_vox = int(lesion.sum())
    les_mm3 = les_vox * vox_mm3

    # ── Per-label computation ────────────────────────────────────────────────────
    labels = np.unique(parc[parc > 0])
    rows = []
    for lbl in labels:
        lbl = int(lbl)
        mask     = parc == lbl
        ovl      = mask & lesion
        lbl_vox  = int(mask.sum())
        ovl_vox  = int(ovl.sum())
        if ovl_vox == 0:
            continue
        lbl_mm3  = lbl_vox * vox_mm3
        ovl_mm3  = ovl_vox * vox_mm3
        pct_lbl  = 100.0 * ovl_vox / lbl_vox  if lbl_vox  > 0 else 0.0
        pct_les  = 100.0 * ovl_vox / les_vox   if les_vox  > 0 else 0.0
        name     = lut.get(lbl, str(lbl))
        rows.append((name, lbl_vox, lbl_mm3, ovl_vox, ovl_mm3, pct_lbl, pct_les))

    rows.sort(key=lambda r: r[6], reverse=True)

    # ── Write text report ────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, 'w') as fh:
        fh.write(f"# KUL_VBG lesion overlap report\n")
        if args.name:
            fh.write(f"# Parcellation: {args.name}\n")
        fh.write(f"# Parc image:   {args.parc}\n")
        fh.write(f"# Lesion mask:  {args.lesion}\n")
        fh.write(f"# LUT:          {args.lut or 'none (label IDs shown)'}\n")
        fh.write(f"# Total lesion: {les_vox} vox | {les_mm3:.1f} mm3\n")
        fh.write(f"# Rows with zero overlap omitted. Sorted by %_of_lesion.\n")
        fh.write("#\n")

        hdr = f"{'Structure':<45} {'Lbl_vox':>8} {'Lbl_mm3':>10} {'Ovl_vox':>8} {'Ovl_mm3':>10} {'%_lbl':>7} {'%_les':>7}\n"
        sep = "-" * len(hdr.rstrip()) + "\n"
        fh.write(hdr)
        fh.write(sep)

        for name, lbl_vox, lbl_mm3, ovl_vox, ovl_mm3, pct_lbl, pct_les in rows:
            fh.write(
                f"{name:<45} {lbl_vox:>8d} {lbl_mm3:>10.1f} "
                f"{ovl_vox:>8d} {ovl_mm3:>10.1f} {pct_lbl:>7.1f} {pct_les:>7.1f}\n"
            )

        fh.write(sep)
        fh.write(f"# {len(rows)} structure(s) with non-zero lesion overlap.\n")

    print(f"Lesion overlap report written: {args.out}  ({len(rows)} overlapping structures)")

    # ── Write / append HTML report ───────────────────────────────────────────────
    if args.html:
        is_first = not os.path.exists(args.html)
        section = _html_section(
            name=args.name or os.path.basename(args.parc),
            rows=rows,
            les_vox=les_vox,
            les_mm3=les_mm3,
            parc_path=args.parc,
            lut_path=args.lut,
            open_tag=is_first,
        )
        append_html_section(args.html, section,
                            les_vox=les_vox, les_mm3=les_mm3,
                            subject=args.subject)
        print(f"HTML report updated: {args.html}")


if __name__ == "__main__":
    main()
