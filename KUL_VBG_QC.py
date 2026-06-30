#!/usr/bin/env python3
"""
KUL_VBG_QC.py — Quality control for the KUL_VBG inpainting pipeline.

Generates:
  - Multi-stage overview PNGs (original → punched → stitched → filled, in MNI + native)
  - Intensity gradient magnitude maps per stage
  - Cross-stage similarity metrics (NCC, SSIM proxy, intensity stats)
  - Artifact report: dark voxels, boundary sharpness, intensity mismatch, stitch score

Usage:
    KUL_VBG_QC.py -s SUBJECT -p PROC_DIR [-o OUTPUT_DIR] [-q QC_DIR] [-m MNI_TEMPLATE]

Arguments:
    -s  Subject ID with sub- prefix (e.g. sub-VBGTest1)
    -p  Full path to proc_VBG/SUBJECT/ directory
    -o  Output directory (for final filled image lookup; defaults to proc dir parent/output_VBG/SUBJECT)
    -q  QC PNG output directory (default: QC/ next to proc dir)
    -m  MNI T1 brain template (auto-detected from FSL/FreeSurfer if not set)
"""

import argparse
import base64
import datetime
import os
import sys
import json
from pathlib import Path

import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.ndimage import uniform_filter, binary_dilation, binary_erosion, zoom


# ── helpers ──────────────────────────────────────────────────────────────────

def load_can(path):
    """Load NIfTI, reorient to RAS canonical, return (data_float32, img)."""
    img = nib.as_closest_canonical(nib.load(str(path)))
    return img.get_fdata(dtype=np.float32), img


def resample_to(source_vol, source_affine, target_shape, target_affine):
    """Trilinear resample source into target voxel grid (numpy only)."""
    tgt = np.zeros(target_shape, dtype=np.float32)
    inv = np.linalg.inv(source_affine)
    for iz in range(target_shape[2]):
        for iy in range(target_shape[1]):
            xs = np.arange(target_shape[0])
            world = target_affine @ np.vstack([
                xs, np.full_like(xs, iy), np.full_like(xs, iz), np.ones_like(xs)])
            vox = (inv @ world)[:3].T
            vx = np.clip(vox[:, 0], 0, source_vol.shape[0] - 1).astype(int)
            vy = np.clip(vox[:, 1], 0, source_vol.shape[1] - 1).astype(int)
            vz = np.clip(vox[:, 2], 0, source_vol.shape[2] - 1).astype(int)
            tgt[:, iy, iz] = source_vol[vx, vy, vz]
    return tgt


def gradient_magnitude(vol, voxel_size=(1.0, 1.0, 1.0)):
    """Voxel-wise gradient magnitude via numpy.gradient."""
    gz, gy, gx = np.gradient(vol.astype(np.float64), *voxel_size)
    return np.sqrt(gx**2 + gy**2 + gz**2).astype(np.float32)


def ncc(a, b, mask=None):
    """Normalized cross-correlation (Pearson r) over mask."""
    if mask is not None:
        a, b = a[mask].astype(np.float64), b[mask].astype(np.float64)
    else:
        a, b = a.ravel().astype(np.float64), b.ravel().astype(np.float64)
    a -= a.mean(); b -= b.mean()
    denom = np.sqrt((a**2).sum() * (b**2).sum())
    return float(np.dot(a, b) / denom) if denom > 1e-12 else 0.0


def ssim_simple(a, b, mask=None):
    """Simplified single-window SSIM over masked region."""
    C1, C2 = (0.01 * 255)**2, (0.03 * 255)**2
    if mask is not None:
        a, b = a[mask].astype(np.float64), b[mask].astype(np.float64)
    else:
        a, b = a.ravel().astype(np.float64), b.ravel().astype(np.float64)
    # normalise to [0, 255] range for C1/C2 constants to make sense
    hi = max(a.max(), b.max(), 1e-6)
    a, b = a / hi * 255, b / hi * 255
    mu_a, mu_b = a.mean(), b.mean()
    sa, sb = a.std(), b.std()
    sab = float(np.cov(a, b)[0, 1]) if len(a) > 1 else 0.0
    num = (2 * mu_a * mu_b + C1) * (2 * sab + C2)
    den = (mu_a**2 + mu_b**2 + C1) * (sa**2 + sb**2 + C2)
    return float(num / den) if abs(den) > 1e-12 else 0.0


def percentile_nz(vol, p, mask=None):
    """Percentile ignoring zeros/near-zeros."""
    v = vol[mask] if mask is not None else vol.ravel()
    v = v[v > 0.001]
    return float(np.percentile(v, p)) if len(v) > 0 else 0.0


def mean_nz(vol, mask=None):
    v = vol[mask] if mask is not None else vol.ravel()
    v = v[v > 0.001]
    return float(v.mean()) if len(v) > 0 else 0.0


def std_nz(vol, mask=None):
    v = vol[mask] if mask is not None else vol.ravel()
    v = v[v > 0.001]
    return float(v.std()) if len(v) > 0 else 0.0


def norm_display(vol, mask=None, plo=2, phi=98):
    """Normalise to [0,1] for display using percentile window from mask/whole vol."""
    ref = vol[mask] if mask is not None else vol.ravel()
    ref = ref[ref > 0.001]
    lo = np.percentile(ref, plo) if len(ref) > 0 else 0.0
    hi = np.percentile(ref, phi) if len(ref) > 0 else 1.0
    return np.clip((vol - lo) / max(hi - lo, 1e-6), 0, 1)


def lesion_com(lmask):
    coords = np.argwhere(lmask)
    if len(coords) == 0:
        return tuple(s // 2 for s in lmask.shape)
    return tuple(coords.mean(axis=0).astype(int))


def find_file(proc_dir, subject, *suffixes):
    """Return first existing file matching any suffix pattern, else None."""
    for suf in suffixes:
        p = proc_dir / f"{subject}_{suf}"
        if p.exists():
            return p
    return None


def find_mni_template():
    candidates = [
        Path(os.environ.get("FSLDIR", "/usr/local/fsl")) / "data/standard/MNI152_T1_1mm_brain.nii.gz",
        Path("/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"),
        Path("/opt/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"),
        Path(os.environ.get("FREESURFER_HOME", "")) / "subjects/cvs_avg35_inMNI152/mri/T1.mgz",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


# ── QC panel rendering ────────────────────────────────────────────────────────

def render_panel(images, lmask, title, out_path, voxel_size=(1,1,1),
                 show_grad=False, show_contour=True):
    """
    images: list of (label, vol_3d) tuples
    lmask: binary lesion mask (may be on a different grid — resampled per-volume)
    Renders axial / coronal / sagittal strips through lesion CoM.
    """
    # Use the first valid volume's shape as the reference for CoM and brain_mask
    ref_vol = next((v for _, v in images if v is not None), None)
    if ref_vol is None:
        print(f"  SKIP (no data): {title}")
        return
    lmask_ref = _match_mask(lmask, ref_vol)
    cx, cy, cz = lesion_com(lmask_ref)
    brain_mask = ref_vol > 0.001

    n = len(images)
    ncols = 3  # ax, cor, sag
    nrows = n

    fig = plt.figure(figsize=(ncols * 3.5, nrows * 3.0), facecolor="black")
    fig.suptitle(title, color="white", fontsize=11, y=1.01)
    gs = gridspec.GridSpec(nrows, ncols, hspace=0.05, wspace=0.05)

    def _slice(vol, ax, plane):
        if plane == "axial":
            s = vol[:, :, cz]
        elif plane == "coronal":
            s = vol[:, cy, :]
        else:
            s = vol[cx, :, :]
        return s.T

    def _contour(ax_obj, mask_sl):
        if not mask_sl.any():
            return
        ax_obj.contour(mask_sl.astype(np.float32), levels=[0.5],
                       colors=["lime"], linewidths=2.5, alpha=0.9)

    planes = ["axial", "coronal", "sagittal"]
    plane_labels = [f"Axial z={cz}", f"Coronal y={cy}", f"Sagittal x={cx}"]

    for row, (label, vol) in enumerate(images):
        if vol is None:
            for col in range(ncols):
                ax = fig.add_subplot(gs[row, col])
                ax.set_facecolor("black")
                ax.text(0.5, 0.5, "not found", color="gray", ha="center", va="center",
                        transform=ax.transAxes, fontsize=7)
                ax.axis("off")
                if row == 0:
                    ax.set_title(plane_labels[col], color="white", fontsize=7)
            fig.text(0.01, (nrows - row - 0.5) / nrows, label, color="yellow",
                     fontsize=7, va="center", rotation=0)
            continue

        lmask_v = _match_mask(lmask, vol)
        bm_v = vol > 0.001
        disp_vol = gradient_magnitude(vol, voxel_size) if show_grad else vol
        disp_vol = norm_display(disp_vol, bm_v)
        cmap = "hot" if show_grad else "gray"

        for col, plane in enumerate(planes):
            ax = fig.add_subplot(gs[row, col])
            sl = _slice(disp_vol, ax, plane)
            ax.imshow(sl, origin="lower", cmap=cmap, vmin=0, vmax=1, aspect="equal")
            if show_contour:
                lsl = _slice(lmask_v.astype(np.uint8), ax, plane)
                _contour(ax, lsl)
            ax.axis("off")
            if row == 0:
                ax.set_title(plane_labels[col], color="white", fontsize=7)

        fig.text(0.01, (nrows - row - 0.5) / nrows, label, color="yellow",
                 fontsize=7, va="center", rotation=0)

    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="black")
    plt.close(fig)
    print(f"  saved: {out_path}")


def render_artifact_panel(vol_filled, lmask, boundary_ring, dark_mask, out_path,
                          voxel_size=(1,1,1)):
    """Three-column artifact map: dark voxels | boundary gradient | intensity diff from native."""
    lmask       = _match_mask(lmask,        vol_filled)
    boundary_ring = _match_mask(boundary_ring, vol_filled) if boundary_ring is not None else lmask
    dark_mask   = _match_mask(dark_mask,    vol_filled) if dark_mask is not None \
                  else np.zeros_like(lmask)
    cx, cy, cz = lesion_com(lmask)
    grad = gradient_magnitude(vol_filled, voxel_size)
    brain_mask = vol_filled > 0.001

    fig, axes = plt.subplots(3, 3, figsize=(10, 10), facecolor="black")
    fig.suptitle("Artifact maps — axial / coronal / sagittal", color="white", fontsize=10)

    planes = ["axial", "coronal", "sagittal"]
    rows_data = [
        ("Dark voxels\n(fill < native P10 WM)", dark_mask.astype(float), "Reds"),
        ("Boundary gradient\nmagnitude", grad * boundary_ring, "hot"),
        ("Fill gradient\n(full lesion)", grad * lmask.astype(float), "hot"),
    ]

    for row, (label, data, cmap) in enumerate(rows_data):
        hi = np.percentile(data[data > 0], 99) if (data > 0).any() else 1.0
        for col, plane in enumerate(planes):
            ax = axes[row, col]
            if plane == "axial":
                sl = data[:, :, cz].T
                bg = norm_display(vol_filled)[:, :, cz].T
            elif plane == "coronal":
                sl = data[:, cy, :].T
                bg = norm_display(vol_filled)[:, cy, :].T
            else:
                sl = data[cx, :, :].T
                bg = norm_display(vol_filled)[cx, :, :].T
            ax.imshow(bg, origin="lower", cmap="gray", vmin=0, vmax=1, aspect="equal")
            if sl.max() > 0:
                ax.imshow(np.clip(sl / max(hi, 1e-6), 0, 1), origin="lower",
                          cmap=cmap, vmin=0, vmax=1, aspect="equal", alpha=0.7)
            ax.axis("off")
            if row == 0:
                ax.set_title(["Axial", "Coronal", "Sagittal"][col],
                             color="white", fontsize=8)
        axes[row, 0].set_ylabel(label, color="yellow", fontsize=7)
        fig.text(0.005, (3 - row - 0.5) / 3, label, color="yellow",
                 fontsize=6, va="center", rotation=90)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="black")
    plt.close(fig)
    print(f"  saved: {out_path}")


def render_metrics_summary(metrics, out_path):
    """Bar charts of per-stage P95, NCC heatmap, artifact scores."""
    stages = [k for k in metrics["per_stage"] if metrics["per_stage"][k]["p95_in_lesion"] > 0]
    p95_vals = [metrics["per_stage"][s]["p95_in_lesion"] for s in stages]
    grad_vals = [metrics["per_stage"][s]["grad_mean_in_lesion"] for s in stages]

    fig = plt.figure(figsize=(14, 8), facecolor="#111111")
    fig.suptitle("KUL_VBG QC — Metric Summary", color="white", fontsize=12)

    # P95 bar
    ax1 = fig.add_subplot(2, 3, 1)
    bars = ax1.barh(stages, p95_vals, color="#4a9eff")
    ax1.set_title("P95 intensity in lesion mask", color="white", fontsize=8)
    ax1.tick_params(colors="white", labelsize=6)
    ax1.set_facecolor("#1a1a1a")
    for spine in ax1.spines.values():
        spine.set_edgecolor("#444")

    # Gradient bar
    ax2 = fig.add_subplot(2, 3, 2)
    ax2.barh(stages, grad_vals, color="#ff7043")
    ax2.set_title("Mean gradient magnitude in lesion", color="white", fontsize=8)
    ax2.tick_params(colors="white", labelsize=6)
    ax2.set_facecolor("#1a1a1a")
    for spine in ax2.spines.values():
        spine.set_edgecolor("#444")

    # NCC heatmap
    sim = metrics.get("similarity", {})
    if sim:
        ax3 = fig.add_subplot(2, 3, 3)
        keys = sorted(sim.keys())
        mat = np.zeros((len(keys), len(keys)))
        labels = []
        for i, k in enumerate(keys):
            parts = k.split("_vs_")
            if len(parts) == 2:
                labels.append(parts[0][:12])
                for j, k2 in enumerate(keys):
                    if sim[k].get("ncc") is not None:
                        mat[i, j] = sim[k]["ncc"]
        im = ax3.imshow(mat, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
        ax3.set_title("NCC (lesion region)", color="white", fontsize=8)
        ax3.tick_params(colors="white", labelsize=5)
        ax3.set_facecolor("#1a1a1a")
        plt.colorbar(im, ax=ax3)

    # Artifact scores
    art = metrics.get("artifacts", {})
    ax4 = fig.add_subplot(2, 1, 2)
    art_labels = []
    art_vals = []
    art_colors = []
    thresholds = {
        "dark_voxel_fraction": 0.05,
        "boundary_grad_ratio": 2.0,
        "wm_intensity_mismatch": 0.10,
        "fill_grad_vs_native_ratio": 1.5,
    }
    for k, v in art.items():
        if isinstance(v, (int, float)):
            art_labels.append(k.replace("_", "\n"))
            art_vals.append(v)
            thr = thresholds.get(k, None)
            art_colors.append("#e53935" if thr is not None and v > thr else "#43a047")
    if art_labels:
        ax4.bar(art_labels, art_vals, color=art_colors)
        for k, thr in thresholds.items():
            if k in art:
                idx = list(art.keys()).index(k)
                if idx < len(art_labels):
                    ax4.axhline(thr, color="orange", linestyle="--", linewidth=0.8,
                                alpha=0.7, label=f"{k} threshold={thr}")
        ax4.set_title("Artifact scores (red = exceeds threshold)", color="white", fontsize=8)
        ax4.tick_params(colors="white", labelsize=6)
        ax4.set_facecolor("#1a1a1a")
        for spine in ax4.spines.values():
            spine.set_edgecolor("#444")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="#111111")
    plt.close(fig)
    print(f"  saved: {out_path}")


# ── metrics ───────────────────────────────────────────────────────────────────

def _match_mask(mask, vol, vol_img=None):
    """Resample mask to vol's shape if they differ. Returns bool array."""
    if mask is None or vol is None:
        return mask
    if mask.shape == vol.shape:
        return mask.astype(bool)
    # resample via nibabel if vol_img provided, else nearest-neighbour zoom
    if vol_img is not None:
        try:
            mask_img = nib.Nifti1Image(mask.astype(np.float32),
                                       vol_img.affine, vol_img.header)
            resampled = resample_from_to(mask_img, vol_img, order=0)
            return np.asarray(resampled.dataobj, dtype=np.float32) > 0.5
        except Exception:
            pass
    # fallback: zoom
    factors = tuple(v / m for v, m in zip(vol.shape, mask.shape))
    return zoom(mask.astype(np.float32), factors, order=0) > 0.5


def compute_per_stage_metrics(label, vol, lmask, native_brain_mask, voxel_size, vol_img=None):
    """Stats within lesion mask and in native brain for a single image."""
    m = {"label": label}
    if vol is None:
        return {k: None for k in ["p95_in_lesion", "mean_in_lesion", "std_in_lesion",
                                   "p95_native", "grad_mean_in_lesion", "grad_std_in_lesion",
                                   "cv_in_lesion"]}
    lbin = _match_mask(lmask, vol, vol_img)
    nbin_raw = _match_mask(native_brain_mask, vol, vol_img) if native_brain_mask is not None else None
    nbin = (nbin_raw & ~lbin) if nbin_raw is not None else ~lbin

    m["p95_in_lesion"]      = percentile_nz(vol, 95, lbin)
    m["mean_in_lesion"]     = mean_nz(vol, lbin)
    m["std_in_lesion"]      = std_nz(vol, lbin)
    m["p95_native"]         = percentile_nz(vol, 95, nbin)
    m["cv_in_lesion"]       = (m["std_in_lesion"] / m["mean_in_lesion"]
                                if m["mean_in_lesion"] > 0 else 0.0)
    grad = gradient_magnitude(vol, voxel_size)
    m["grad_mean_in_lesion"] = mean_nz(grad, lbin)
    m["grad_std_in_lesion"]  = std_nz(grad, lbin)
    m["grad_mean_native"]    = mean_nz(grad, nbin)
    return m


def compute_similarity(label_a, vol_a, label_b, vol_b, lmask):
    """NCC and SSIM proxy between two volumes inside the lesion mask."""
    key = f"{label_a}_vs_{label_b}"
    if vol_a is None or vol_b is None:
        return key, {"ncc": None, "ssim": None}
    lbin = lmask.astype(bool)
    # ensure matching shape — skip if mismatch (different spaces)
    if vol_a.shape != vol_b.shape:
        return key, {"ncc": None, "ssim": None, "note": "shape mismatch"}
    return key, {
        "ncc":  ncc(vol_a, vol_b, lbin),
        "ssim": ssim_simple(vol_a, vol_b, lbin),
    }


def compute_artifacts(vol_filled, vol_native_clean, lmask, boundary_ring,
                      voxel_size=(1,1,1)):
    """
    Returns artifact metrics dict:
      dark_voxel_fraction      — fraction of lesion voxels below native WM P10
      wm_intensity_mismatch    — |P95(fill) - P95(native)| / P95(native)
      boundary_grad_ratio      — mean gradient at boundary ring / mean gradient native
      fill_grad_vs_native_ratio— mean gradient inside lesion / mean gradient native
      stitch_score             — 99th percentile of gradient at boundary ring (normalised)
      n_dark_voxels            — absolute count
    """
    art = {}
    lbin = lmask.astype(bool)
    nbin = ~lbin & (vol_native_clean > 0.001) if vol_native_clean is not None else ~lbin

    # Dark voxels: fill intensity < native P10 WM (P10 of native outside lesion)
    if vol_native_clean is not None:
        native_p10 = percentile_nz(vol_native_clean, 10, nbin)
        native_p95 = percentile_nz(vol_native_clean, 95, nbin)
        fill_p95   = percentile_nz(vol_filled, 95, lbin)
        dark_mask  = lbin & (vol_filled < native_p10)
        art["n_dark_voxels"]        = int(dark_mask.sum())
        art["dark_voxel_fraction"]  = float(dark_mask.sum() / max(lbin.sum(), 1))
        art["wm_intensity_mismatch"]= float(abs(fill_p95 - native_p95) / max(native_p95, 1e-6))
    else:
        dark_mask = np.zeros_like(lbin)
        native_p10 = 0.0
        art["n_dark_voxels"]         = 0
        art["dark_voxel_fraction"]   = 0.0
        art["wm_intensity_mismatch"] = 0.0

    grad = gradient_magnitude(vol_filled, voxel_size)
    ring_bin = boundary_ring.astype(bool) if boundary_ring is not None else lbin
    native_grad_mean = mean_nz(grad, nbin)

    art["boundary_grad_ratio"]      = (mean_nz(grad, ring_bin) / max(native_grad_mean, 1e-6))
    art["fill_grad_vs_native_ratio"]= (mean_nz(grad, lbin) / max(native_grad_mean, 1e-6))
    art["stitch_score"]             = (percentile_nz(grad, 99, ring_bin) /
                                       max(percentile_nz(grad, 99, nbin), 1e-6))

    return art, dark_mask


# ── HTML report ──────────────────────────────────────────────────────────────

_QC_CSS = """
body{font-family:Arial,sans-serif;font-size:13px;margin:20px;background:#f8f9fa;color:#222}
h1{color:#2c3e50;border-bottom:2px solid #3498db;padding-bottom:6px}
h2{color:#34495e;font-size:14px;margin:18px 0 6px}
.meta{background:#fff;border:1px solid #dee2e6;border-radius:4px;padding:10px 14px;
      margin-bottom:16px;font-size:12px;color:#555}
details{background:#fff;border:1px solid #dee2e6;border-radius:4px;margin-bottom:10px}
summary{cursor:pointer;padding:10px 14px;font-weight:bold;font-size:13px;
        background:#f1f3f5;border-radius:4px;list-style:none;user-select:none}
summary::before{content:"▶ ";font-size:10px;margin-right:4px}
details[open] summary::before{content:"▼ "}
.section-body{padding:10px 14px 14px}
.checklist{list-style:none;padding:0;margin:0}
.checklist li{padding:4px 0;font-size:12px;border-bottom:1px solid #f0f0f0}
.badge{display:inline-block;border-radius:3px;padding:1px 7px;font-size:11px;
       font-weight:bold;color:#fff;margin-right:6px;min-width:40px;text-align:center}
.pass{background:#27ae60}.warn{background:#f39c12}.fail{background:#e74c3c}.skip{background:#95a5a6}
table{border-collapse:collapse;width:100%;font-size:12px;margin-top:6px}
th{background:#3498db;color:#fff;padding:5px 8px;text-align:left}
td{padding:4px 8px;border-bottom:1px solid #eee}
tr:hover td{background:#f0f7ff}
.qc-img{max-width:100%;border:1px solid #dee2e6;border-radius:4px;margin:6px 0;display:block}
.artifact-ok{color:#27ae60;font-weight:bold}
.artifact-fail{color:#e74c3c;font-weight:bold}
.overall-pass{background:#d4edda;border:1px solid #c3e6cb;color:#155724;
              border-radius:4px;padding:10px 14px;margin-top:10px;font-weight:bold}
.overall-fail{background:#f8d7da;border:1px solid #f5c6cb;color:#721c24;
              border-radius:4px;padding:10px 14px;margin-top:10px;font-weight:bold}
a{color:#3498db}
"""

_THRESHOLDS = {
    "dark_voxel_fraction": 0.05,
    "boundary_grad_ratio": 2.0,
    "wm_intensity_mismatch": 0.10,
    "fill_grad_vs_native_ratio": 1.5,
    "stitch_score": 2.5,
}


def _png_b64(path):
    if path is None or not Path(path).exists():
        return None
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


def _badge(status):
    cls = {"PASS": "pass", "WARN": "warn", "FAIL": "fail", "SKIP": "skip"}.get(status, "skip")
    return f'<span class="badge {cls}">{status}</span>'


def write_qc_html(html_path, subject, proc, out_d, scripts_dir,
                  metrics, sim_pairs, art_metrics, qc_dir, overlap_html=None):
    date_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    # ── checklist ────────────────────────────────────────────────────────────
    def done(fname, parent=None):
        p = Path(parent or scripts_dir) / fname if scripts_dir else None
        if p is None:
            return "SKIP"
        return "PASS" if p.exists() else "FAIL"

    def exists(path):
        return "PASS" if (path is not None and Path(path).exists()) else "FAIL"

    def non_empty(path):
        if path is None or not Path(path).exists():
            return "FAIL"
        return "PASS" if Path(path).stat().st_size > 0 else "WARN"

    lmask_path = next((proc / f"{subject}_{s}"
                       for s in ("Lmask_in_T1_bin.nii.gz", "Lmask_in_T1_bins3.nii.gz")
                       if (proc / f"{subject}_{s}").exists()), None)
    fill_path  = next((out_d / f"{subject}_{s}"
                       for s in ("T1_stdOri_filled_1.nii.gz",)
                       if (out_d / f"{subject}_{s}").exists()), None)
    brain_path = next((proc / f"{subject}_{s}"
                       for s in ("Brain_clean.nii.gz", "T1_brain_clean.nii.gz")
                       if (proc / f"{subject}_{s}").exists()), None)

    checklist = [
        ("Lesion mask",              non_empty(lmask_path)),
        ("Brain mask / clean T1",    non_empty(brain_path)),
        ("VBG-filled image",         non_empty(fill_path)),
        ("recon-all",                done("recon-all.done")),
        ("Lausanne/Glasser parc",    done("multiscale_parc.done")),
        ("Thalamic nuclei",          done("thalamic_nuclei.done")),
        ("Brainstem substructures",  done("brainstem_subregions.done")),
        ("Hippo/amygdala subregions",done("hippo_amygdala.done")),
        ("Hypothalamic subunits",    done("hypothalamic_subunits.done")),
    ]

    check_html = "<ul class='checklist'>"
    for label, status in checklist:
        check_html += f"<li>{_badge(status)} {label}</li>"
    check_html += "</ul>"

    # ── images ───────────────────────────────────────────────────────────────
    png_labels = [
        (f"{subject}_QC_01_mni_stages.png",      "MNI-space pipeline stages"),
        (f"{subject}_QC_02_mni_gradients.png",   "MNI-space gradient maps"),
        (f"{subject}_QC_03_native_stages.png",   "Native-space fill stages"),
        (f"{subject}_QC_04_native_gradients.png","Native-space gradient maps"),
        (f"{subject}_QC_05_artifacts.png",       "Artifact maps"),
        (f"{subject}_QC_06_metrics_summary.png", "Metrics summary"),
    ]
    img_html = ""
    for fname, label in png_labels:
        b64 = _png_b64(qc_dir / fname)
        if b64:
            img_html += (f"<p style='margin:8px 0 2px;font-weight:bold;font-size:12px'>{label}</p>"
                         f"<img class='qc-img' src='data:image/png;base64,{b64}' alt='{label}'>")
        else:
            img_html += f"<p style='color:#999;font-size:12px'>{label} — not generated</p>"

    # ── metrics table ─────────────────────────────────────────────────────────
    stage_rows = ""
    for label, m in metrics.get("per_stage", {}).items():
        if m and m.get("p95_in_lesion") is not None:
            stage_rows += (
                f"<tr><td>{label}</td>"
                f"<td>{m['p95_in_lesion']:.2f}</td>"
                f"<td>{m.get('mean_in_lesion', 0):.2f}</td>"
                f"<td>{m.get('grad_mean_in_lesion', 0):.3f}</td>"
                f"<td>{m.get('cv_in_lesion', 0):.3f}</td></tr>"
            )
    metrics_html = f"""
<table>
<thead><tr><th>Stage</th><th>P95 in lesion</th><th>Mean in lesion</th>
<th>Grad mean</th><th>CV</th></tr></thead>
<tbody>{stage_rows}</tbody>
</table>"""

    # ── similarity table ──────────────────────────────────────────────────────
    sim_rows = ""
    for key, res in sim_pairs:
        ncc_v  = f"{res['ncc']:.4f}"  if res.get("ncc")  is not None else "N/A"
        ssim_v = f"{res['ssim']:.4f}" if res.get("ssim") is not None else "N/A"
        sim_rows += f"<tr><td>{key}</td><td>{ncc_v}</td><td>{ssim_v}</td></tr>"
    sim_html = (f"<table><thead><tr><th>Pair</th><th>NCC</th><th>SSIM</th></tr></thead>"
                f"<tbody>{sim_rows}</tbody></table>") if sim_rows else "<p>No pairs computed.</p>"

    # ── artifact scores ────────────────────────────────────────────────────────
    art_rows = ""
    fails = []
    for k, v in art_metrics.items():
        if isinstance(v, float):
            thr = _THRESHOLDS.get(k)
            over = thr is not None and v > thr
            cls = "artifact-fail" if over else "artifact-ok"
            flag = f" &nbsp;<b>↑ &gt; {thr}</b>" if over else ""
            art_rows += f"<tr><td>{k}</td><td class='{cls}'>{v:.4f}{flag}</td></tr>"
            if over:
                fails.append(f"{k} = {v:.3f} (threshold {thr})")
        else:
            art_rows += f"<tr><td>{k}</td><td>{v}</td></tr>"
    art_html = (f"<table><thead><tr><th>Metric</th><th>Value</th></tr></thead>"
                f"<tbody>{art_rows}</tbody></table>") if art_rows else "<p>Not computed.</p>"

    # ── overall verdict ────────────────────────────────────────────────────────
    checklist_fails = [label for label, s in checklist if s == "FAIL"]
    all_fails = fails + [f"Pipeline step missing: {l}" for l in checklist_fails]
    if all_fails:
        verdict = ("<div class='overall-fail'>&#10008; FAIL<ul>"
                   + "".join(f"<li>{f}</li>" for f in all_fails)
                   + "</ul></div>")
    else:
        verdict = "<div class='overall-pass'>&#10004; PASS — no thresholds exceeded, all steps done</div>"

    # ── overlap report link ────────────────────────────────────────────────────
    overlap_section = ""
    if overlap_html and Path(overlap_html).exists():
        overlap_section = (
            f"<details><summary>Lesion Overlap Report</summary>"
            f"<div class='section-body'>"
            f"<p><a href='{overlap_html}' target='_blank'>Open lesion overlap report</a></p>"
            f"</div></details>"
        )

    bilateral = metrics.get("bilateral", False)
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>KUL_VBG QC — {subject}</title>
<style>{_QC_CSS}</style>
</head>
<body>
<h1>KUL_VBG QC Report</h1>
<div class="meta">
  <b>Subject:</b> {subject} &nbsp;|&nbsp;
  <b>Date:</b> {date_str} &nbsp;|&nbsp;
  <b>Laterality:</b> {'Bilateral' if bilateral else 'Unilateral'} &nbsp;|&nbsp;
  <b>Proc dir:</b> {proc}
</div>

{verdict}

<details open>
<summary>Pipeline checklist</summary>
<div class="section-body">{check_html}</div>
</details>

<details open>
<summary>QC images</summary>
<div class="section-body">{img_html}</div>
</details>

<details>
<summary>Per-stage intensity metrics</summary>
<div class="section-body">{metrics_html}</div>
</details>

<details>
<summary>Cross-stage similarity (NCC / SSIM in lesion)</summary>
<div class="section-body">{sim_html}</div>
</details>

<details>
<summary>Artifact detection</summary>
<div class="section-body">{art_html}</div>
</details>

{overlap_section}

</body>
</html>
"""
    _hp = Path(html_path)
    _hp.parent.mkdir(parents=True, exist_ok=True)
    _hp.write_text(html)
    print(f"  QC HTML report: {_hp}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="KUL_VBG QC script")
    parser.add_argument("-s", "--subject", required=True,
                        help="Subject ID with sub- prefix (e.g. sub-VBGTest1)")
    parser.add_argument("-p", "--proc_dir", required=True,
                        help="Path to proc_VBG/SUBJECT/ directory")
    parser.add_argument("-o", "--output_dir", default=None,
                        help="Path to output_VBG/SUBJECT/ directory")
    parser.add_argument("-q", "--qc_dir", default=None,
                        help="QC PNG output directory (default: QC/ next to proc_dir)")
    parser.add_argument("-m", "--mni_template", default=None,
                        help="MNI T1 brain template (auto-detected if not set)")
    parser.add_argument("-d", "--scripts_dir", default=None,
                        help="FreeSurfer subject scripts/ dir where *.done files live")
    parser.add_argument("--html", default=None,
                        help="Write self-contained QC HTML report to this path")
    parser.add_argument("--overlap_html", default=None,
                        help="Path to lesion overlap HTML report (linked from QC report)")
    args = parser.parse_args()

    S_raw = args.subject  # as given (may or may not have sub- prefix)
    proc_base = Path(args.proc_dir)
    if not proc_base.exists():
        sys.exit(f"ERROR: proc dir not found: {proc_base}")

    # ── auto-detect subject subdir and file prefix ────────────────────────────
    # KUL_VBG writes files as {subj_dir}/{subj_prefix}_{suffix}.
    # The user may pass the parent (proc_VBG/) or the subdir (proc_VBG/sub-X/).
    # Files may be prefixed sub-{S} or just {S}.
    def _probe(directory, prefix):
        """True if directory exists and contains at least one {prefix}_*.nii.gz."""
        d = Path(directory)
        return d.is_dir() and any(d.glob(f"{prefix}_*.nii.gz"))

    proc = None
    S = None
    for candidate_prefix in [S_raw, f"sub-{S_raw}"]:
        # user gave the subdir directly
        if _probe(proc_base, candidate_prefix):
            proc, S = proc_base, candidate_prefix
            break
        # user gave the parent — look one level deeper
        for subdir in [proc_base / candidate_prefix, proc_base / S_raw, proc_base / f"sub-{S_raw}"]:
            if _probe(subdir, candidate_prefix):
                proc, S = subdir, candidate_prefix
                break
        if proc:
            break

    if proc is None:
        # fallback: use as-given so error messages are informative
        proc, S = proc_base, S_raw
        print(f"  WARN: could not auto-detect files under {proc_base} for subject "
              f"'{S_raw}' — proceeding anyway (expect MISS entries)")

    # ── output dir ────────────────────────────────────────────────────────────
    if args.output_dir:
        out_d_base = Path(args.output_dir)
        # same auto-detect: user may give parent or subdir
        if _probe(out_d_base, S):
            out_d = out_d_base
        elif _probe(out_d_base / S, S):
            out_d = out_d_base / S
        else:
            out_d = out_d_base
    else:
        out_d = proc.parent.parent / "output_VBG" / S
    out_d.mkdir(parents=True, exist_ok=True)

    qc_dir = Path(args.qc_dir) if args.qc_dir else proc.parent / "QC"
    qc_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nKUL_VBG QC  subject={S}")
    print(f"  proc   : {proc}")
    print(f"  output : {out_d}")
    print(f"  qc_dir : {qc_dir}\n")

    # ── file discovery ────────────────────────────────────────────────────────
    def F(*suf):
        return find_file(proc, S, *suf)

    mni_tpl_path = Path(args.mni_template) if args.mni_template else find_mni_template()

    def O(*suf):
        """find_file but searches output dir."""
        return find_file(out_d, S, *suf)

    files = {
        # MNI1-space images (all live in proc dir)
        "t1_mni1":      F("T1_brain_inMNI1_Warped.nii.gz"),
        "punched":      F("T1brain_inMNI1_punched.nii.gz"),
        "punched_norm": F("T1brain_inMNI1_punched_norm.nii.gz"),
        "stitched":     F("tmp_s2T1_CSFGMCBWM.nii.gz"),
        "filled_mni":   F("T1_brain_Temp_bil_Lmask_filled1.nii.gz"),
        "sti2fill":     F("stitchT12filled_brain_Warped.nii.gz"),
        # lesion masks
        "lmask_mni1":   F("Lmask_bin_s3_inMNI1.nii.gz", "Lmask_bin_inMNI1.nii.gz"),
        "lmask_native": F("Lmask_in_T1_bins3.nii.gz", "Lmask_in_T1_bin.nii.gz",
                          "L_mask_in_T1_bin.nii.gz"),
        "boundary_ring":F("lesion_boundary_ring.nii.gz"),
        # native-space images (proc dir first, then output dir)
        "t1_clean":     F("Brain_clean.nii.gz", "T1_brain_clean.nii.gz"),
        "initial_fill": F("T1_brain_bk2anat1_InverseWarped.nii.gz",
                          "T1_initial_filled_brain.nii.gz") or \
                        O("T1_initial_filled_brain.nii.gz"),
        "final_fill":   F("T1_stdOri_filled_1.nii.gz") or \
                        O("T1_stdOri_filled.nii.gz", "T1_nat_filled_brain.nii.gz"),
        # MNI template
        "mni_template": mni_tpl_path,
    }

    print("File discovery:")
    for k, v in files.items():
        status = "OK  " if (v is not None and Path(v).exists()) else "MISS"
        print(f"  [{status}] {k:20s}  {v if v else '—'}")
    print()

    # ── load images ───────────────────────────────────────────────────────────
    def L(key):
        p = files[key]
        if p is None or not Path(p).exists():
            return None, None
        try:
            v, img = load_can(p)
            return v, img
        except Exception as e:
            print(f"  WARN: could not load {key}: {e}")
            return None, None

    t1_mni1,    t1_mni1_img    = L("t1_mni1")
    punched,    punched_img    = L("punched")
    stitched,   stitched_img   = L("stitched")
    filled_mni, filled_mni_img = L("filled_mni")
    sti2fill,   sti2fill_img   = L("sti2fill")
    lmask_mni1, lmask_mni1_img = L("lmask_mni1")
    lmask_nat,  lmask_nat_img  = L("lmask_native")
    bnd_ring,   _              = L("boundary_ring")
    t1_clean,   t1_clean_img   = L("t1_clean")
    init_fill,  init_fill_img  = L("initial_fill")
    final_fill, final_fill_img = L("final_fill")
    mni_tpl,    mni_tpl_img    = L("mni_template")

    # binary masks
    lmask_mni1_bin = (lmask_mni1 > 0.5).astype(bool) if lmask_mni1 is not None else None
    lmask_nat_bin  = (lmask_nat  > 0.5).astype(bool) if lmask_nat  is not None else None
    bnd_ring_bin   = (bnd_ring   > 0.5).astype(bool) if bnd_ring   is not None else None

    if lmask_mni1_bin is None and lmask_nat_bin is None:
        sys.exit("ERROR: no lesion mask found — cannot compute metrics.")

    # voxel sizes
    def vox(img):
        if img is None:
            return (1.0, 1.0, 1.0)
        return tuple(float(v) for v in img.header.get_zooms()[:3])

    vs_mni  = vox(t1_mni1_img)
    vs_nat  = vox(t1_clean_img if t1_clean_img else final_fill_img)

    # native brain mask (non-lesion)
    nat_brain_mask = (t1_clean > 0.001) if t1_clean is not None else None

    # ── per-stage metrics ─────────────────────────────────────────────────────
    print("Computing per-stage metrics…")
    mni_stages = [
        ("T1_inMNI1",    t1_mni1,    lmask_mni1_bin, None,          vs_mni, t1_mni1_img),
        ("Punched",      punched,    lmask_mni1_bin, None,          vs_mni, punched_img),
        ("Stitched",     stitched,   lmask_mni1_bin, None,          vs_mni, stitched_img),
        ("Filled_MNI",   filled_mni, lmask_mni1_bin, None,          vs_mni, filled_mni_img),
        ("Stitch2Fill",  sti2fill,   lmask_mni1_bin, None,          vs_mni, sti2fill_img),
        ("MNI_template", mni_tpl,    lmask_mni1_bin, None,          vs_mni, mni_tpl_img),
    ]
    nat_stages = [
        ("T1_clean",     t1_clean,   lmask_nat_bin,  nat_brain_mask, vs_nat, t1_clean_img),
        ("Initial_fill", init_fill,  lmask_nat_bin,  nat_brain_mask, vs_nat, init_fill_img),
        ("Final_fill",   final_fill, lmask_nat_bin,  nat_brain_mask, vs_nat, final_fill_img),
    ]
    per_stage = {}
    for label, vol, lm, nbm, vs, vimg in mni_stages + nat_stages:
        if lm is None:
            continue
        m = compute_per_stage_metrics(label, vol, lm, nbm, vs, vol_img=vimg)
        per_stage[label] = m
        if vol is not None:
            print(f"  {label:20s}  P95={m['p95_in_lesion']:.2f}  "
                  f"mean={m['mean_in_lesion']:.2f}  "
                  f"grad={m['grad_mean_in_lesion']:.3f}")

    # ── cross-stage similarity ─────────────────────────────────────────────────
    print("\nComputing cross-stage similarity…")
    sim_pairs = []
    # MNI-space pairs — only compare volumes on the same grid
    if lmask_mni1_bin is not None:
        mni_vols = [("T1_inMNI1", t1_mni1, t1_mni1_img),
                    ("Punched",   punched,  punched_img),
                    ("Stitched",  stitched, stitched_img),
                    ("Filled_MNI", filled_mni, filled_mni_img)]
        for i, (la, va, _) in enumerate(mni_vols):
            for lb, vb, vb_img in mni_vols[i+1:]:
                if va is None or vb is None or va.shape != vb.shape:
                    continue
                lm = _match_mask(lmask_mni1_bin, va, t1_mni1_img)
                key, res = compute_similarity(la, va, lb, vb, lm)
                sim_pairs.append((key, res))
                if res["ncc"] is not None:
                    print(f"  {key:40s}  NCC={res['ncc']:.4f}  SSIM={res['ssim']:.4f}")
    # native-space pairs
    if lmask_nat_bin is not None:
        nat_vols = [("T1_clean",     t1_clean,   t1_clean_img),
                    ("Initial_fill", init_fill,  init_fill_img),
                    ("Final_fill",   final_fill, final_fill_img)]
        for i, (la, va, va_img) in enumerate(nat_vols):
            for lb, vb, _ in nat_vols[i+1:]:
                if va is None or vb is None or va.shape != vb.shape:
                    continue
                lm = _match_mask(lmask_nat_bin, va, va_img)
                key, res = compute_similarity(la, va, lb, vb, lm)
                sim_pairs.append((key, res))
                if res["ncc"] is not None:
                    print(f"  {key:40s}  NCC={res['ncc']:.4f}  SSIM={res['ssim']:.4f}")

    # ── artifact detection ─────────────────────────────────────────────────────
    print("\nRunning artifact detection…")
    art_metrics = {}
    dark_mask_nat = None
    if final_fill is not None and lmask_nat_bin is not None:
        _lm_ff = _match_mask(lmask_nat_bin, final_fill, final_fill_img)
        _bnd   = _match_mask(bnd_ring_bin,  final_fill, final_fill_img)
        _tc    = t1_clean if (t1_clean is not None and t1_clean.shape == final_fill.shape) else None
        art_metrics, dark_mask_nat = compute_artifacts(
            final_fill, _tc, _lm_ff, _bnd, vs_nat)
        print(f"  dark_voxel_fraction   = {art_metrics['dark_voxel_fraction']:.4f}"
              + (" ← ABOVE THRESHOLD" if art_metrics["dark_voxel_fraction"] > 0.05 else ""))
        print(f"  wm_intensity_mismatch = {art_metrics['wm_intensity_mismatch']:.4f}"
              + (" ← ABOVE THRESHOLD" if art_metrics["wm_intensity_mismatch"] > 0.10 else ""))
        print(f"  boundary_grad_ratio   = {art_metrics['boundary_grad_ratio']:.4f}"
              + (" ← ABOVE THRESHOLD (possible stitch line)"
                 if art_metrics["boundary_grad_ratio"] > 2.0 else ""))
        print(f"  stitch_score          = {art_metrics['stitch_score']:.4f}"
              + (" ← HIGH (check boundary)" if art_metrics["stitch_score"] > 2.5 else ""))
        print(f"  n_dark_voxels         = {art_metrics['n_dark_voxels']}")
    else:
        print("  SKIP: final fill or native lesion mask not found")

    # ── detect bilateral lesion ────────────────────────────────────────────────
    bilateral = False
    if lmask_mni1_bin is not None:
        mid = lmask_mni1_bin.shape[0] // 2
        total = lmask_mni1_bin.sum()
        side_a = lmask_mni1_bin[:mid].sum()
        side_b = lmask_mni1_bin[mid:].sum()
        # require at least 5% of total lesion on each side to call bilateral
        bilateral = (total > 0 and
                     side_a / total >= 0.05 and
                     side_b / total >= 0.05)
        pct_a = 100.0 * side_a / total if total > 0 else 0
        pct_b = 100.0 * side_b / total if total > 0 else 0
        print(f"\n  Lesion laterality: {'BILATERAL' if bilateral else 'UNILATERAL'} "
              f"(side_a={pct_a:.1f}%  side_b={pct_b:.1f}%)")

    # ── assemble metrics dict ──────────────────────────────────────────────────
    metrics = {
        "subject": S,
        "bilateral": bilateral,
        "per_stage": per_stage,
        "similarity": dict(sim_pairs),
        "artifacts": art_metrics,
    }
    metrics_path = qc_dir / f"{S}_QC_metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2, default=lambda x: None)
    print(f"\n  metrics saved: {metrics_path}")

    # ── render PNGs ───────────────────────────────────────────────────────────
    print("\nRendering QC panels…")

    # Panel 1: MNI-space multi-stage overview
    # Note: mni_tpl is in FSL MNI152 space (182×218×182), not the pipeline's FS MNI1 space
    # (256³) — shown separately to avoid affine mismatch with the lesion overlay.
    if lmask_mni1_bin is not None and t1_mni1 is not None:
        mni_panels = [
            ("T1 original (MNI1)",        t1_mni1),
            ("Punched (lesion=0)",         punched),
            ("Stitched composite",         stitched),
            ("Stitch→fill aligned",        sti2fill),
            ("Filled MNI (rough splice)",  filled_mni),
        ]
        render_panel(mni_panels, lmask_mni1_bin,
                     f"{S} — MNI1-space pipeline stages",
                     qc_dir / f"{S}_QC_01_mni_stages.png",
                     voxel_size=vs_mni, show_grad=False)

    # Panel 1b: MNI reference template (different space — no lesion overlay)
    if mni_tpl is not None and lmask_mni1_bin is not None:
        render_panel([("MNI template (FSL MNI152)", mni_tpl)],
                     lmask_mni1_bin,
                     f"{S} — MNI152 reference template",
                     qc_dir / f"{S}_QC_01b_mni_template.png",
                     voxel_size=vs_mni, show_contour=False)

    # Panel 2: MNI-space gradient maps
    if lmask_mni1_bin is not None and t1_mni1 is not None:
        render_panel(mni_panels, lmask_mni1_bin,
                     f"{S} — MNI1-space gradient magnitude maps",
                     qc_dir / f"{S}_QC_02_mni_gradients.png",
                     voxel_size=vs_mni, show_grad=True)

    # Panel 3: Native-space stages
    if lmask_nat_bin is not None:
        nat_panels = [
            ("T1 clean (lesion=0)", t1_clean),
            ("Initial fill (InvWarp)", init_fill),
            ("Final fill (Task 2.5)", final_fill),
        ]
        render_panel(nat_panels, lmask_nat_bin,
                     f"{S} — Native-space fill stages",
                     qc_dir / f"{S}_QC_03_native_stages.png",
                     voxel_size=vs_nat, show_grad=False)

    # Panel 4: Native-space gradient maps
    if lmask_nat_bin is not None:
        render_panel(nat_panels, lmask_nat_bin,
                     f"{S} — Native-space gradient magnitude maps",
                     qc_dir / f"{S}_QC_04_native_gradients.png",
                     voxel_size=vs_nat, show_grad=True)

    # Panel 5: Artifact maps
    if final_fill is not None and lmask_nat_bin is not None:
        _bnd_disp = bnd_ring_bin if bnd_ring_bin is not None else (
            (binary_dilation(lmask_nat_bin, iterations=2).astype(np.uint8) -
             binary_erosion(lmask_nat_bin, iterations=2).astype(np.uint8)).astype(bool))
        _dark_disp = dark_mask_nat if dark_mask_nat is not None else np.zeros_like(lmask_nat_bin)
        render_artifact_panel(final_fill, lmask_nat_bin, _bnd_disp, _dark_disp,
                              qc_dir / f"{S}_QC_05_artifacts.png",
                              voxel_size=vs_nat)

    # Panel 6: Metrics summary chart
    render_metrics_summary(metrics, qc_dir / f"{S}_QC_06_metrics_summary.png")

    # ── text summary report ────────────────────────────────────────────────────
    report_path = qc_dir / f"{S}_QC_report.txt"
    with open(report_path, "w") as rp:
        def W(line=""):
            print(line); rp.write(line + "\n")

        W("=" * 70)
        W(f"KUL_VBG QC REPORT  —  {S}")
        W("=" * 70)
        W(f"Lesion: {'BILATERAL' if bilateral else 'UNILATERAL'}")
        W()
        W("── Per-stage intensity (P95 in lesion, gradient mean in lesion) ──")
        for label, m in per_stage.items():
            if m.get("p95_in_lesion") is None:
                W(f"  {label:20s}  [not available]")
            else:
                W(f"  {label:20s}  P95={m['p95_in_lesion']:.2f}  "
                  f"P95_native={m.get('p95_native', 0):.2f}  "
                  f"grad={m['grad_mean_in_lesion']:.3f}  "
                  f"CV={m.get('cv_in_lesion', 0):.3f}")
        W()
        W("── Cross-stage similarity (NCC / SSIM in lesion mask) ──")
        for key, res in sim_pairs:
            ncc_v = f"{res['ncc']:.4f}" if res["ncc"] is not None else "N/A"
            ssim_v = f"{res['ssim']:.4f}" if res["ssim"] is not None else "N/A"
            W(f"  {key:45s}  NCC={ncc_v}  SSIM={ssim_v}")
        W()
        W("── Artifact detection ──")
        if art_metrics:
            for k, v in art_metrics.items():
                if isinstance(v, float):
                    flag = ""
                    thresh = {"dark_voxel_fraction": 0.05,
                              "boundary_grad_ratio": 2.0,
                              "wm_intensity_mismatch": 0.10,
                              "fill_grad_vs_native_ratio": 1.5,
                              "stitch_score": 2.5}.get(k)
                    if thresh is not None and v > thresh:
                        flag = "  ← ABOVE THRESHOLD"
                    W(f"  {k:35s} = {v:.4f}{flag}")
                else:
                    W(f"  {k:35s} = {v}")
        else:
            W("  [not computed]")
        W()
        W("── Overall pass/fail ──")
        fails = []
        if art_metrics.get("dark_voxel_fraction", 0) > 0.05:
            fails.append(f"dark_voxel_fraction={art_metrics['dark_voxel_fraction']:.3f} > 0.05")
        if art_metrics.get("wm_intensity_mismatch", 0) > 0.10:
            fails.append(f"wm_intensity_mismatch={art_metrics['wm_intensity_mismatch']:.3f} > 0.10")
        if art_metrics.get("boundary_grad_ratio", 0) > 2.0:
            fails.append(f"boundary_grad_ratio={art_metrics['boundary_grad_ratio']:.3f} > 2.0 (stitch line?)")
        if fails:
            W("  FAIL:")
            for f in fails:
                W(f"    - {f}")
        else:
            W("  PASS — no artifact thresholds exceeded")
        W()
        W(f"Metrics JSON : {metrics_path}")
        W(f"QC PNGs      : {qc_dir}/")
        W("=" * 70)

    print(f"\n  report saved: {report_path}")

    # ── HTML report ────────────────────────────────────────────────────────────
    scripts_dir = Path(args.scripts_dir) if args.scripts_dir else None
    if args.html:
        _hp = Path(args.html)
        if _hp.is_dir():
            html_out = _hp / f"{S}_QC_report.html"
        elif _hp.is_absolute():
            html_out = _hp
        else:
            html_out = qc_dir / _hp.name
    else:
        html_out = qc_dir / f"{S}_QC_report.html"
    # safety net: if resolved path is still a directory, put filename inside it
    if html_out.exists() and html_out.is_dir():
        html_out = html_out / f"{S}_QC_report.html"
    write_qc_html(
        html_path    = html_out,
        subject      = S,
        proc         = proc,
        out_d        = out_d,
        scripts_dir  = scripts_dir,
        metrics      = metrics,
        sim_pairs    = sim_pairs,
        art_metrics  = art_metrics,
        qc_dir       = qc_dir,
        overlap_html = args.overlap_html,
    )

    print("\nDone.\n")


if __name__ == "__main__":
    main()
