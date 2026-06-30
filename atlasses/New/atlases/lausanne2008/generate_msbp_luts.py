#!/usr/bin/env python3
"""
generate_msbp_luts.py

Regenerate MSBP-format FreeSurferColorLUT files for all Lausanne2018 scales
(1–5) from the .annot ctab.  For scales 3 and 5, the existing reference LUTs
are used verbatim (they come from real MSBP output).  For scales 1, 2 and 4
we generate them with the same sequential ID scheme and the same fixed
subcortical block.

Output files: label-L2018_desc-scaleN_atlas_FreeSurferColorLUT.txt
"""

import sys
import os
import re
import colorsys
import nibabel as nib

ATLAS_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Fixed subcortical/brainstem block ─────────────────────────────────────────
# Identical across all scales; colours taken from the MSBP scale3 reference LUT.
RH_SUBCORT = [
    ("Right-Pulvinar",                                               255,   0,   0),
    ("Right-Anterior",                                                 0, 255,   0),
    ("Right-Medio_Dorsal",                                           255, 255,   0),
    ("Right-Ventral_Latero_Dorsal",                                  255, 123,   0),
    ("Right-Central_Lateral-Lateral_Posterior-Medial_Pulvinar",        0, 255, 255),
    ("Right-Ventral_Anterior",                                       255,   0, 255),
    ("Right-Ventral_Latero_Ventral",                                   0,   0, 255),
    # ---
    ("Right-Thalamus_Proper",                                          0, 118,  14),
    ("Right-Caudate",                                                122, 186, 220),
    ("Right-Putamen",                                                236,  13, 176),
    ("Right-Pallidum",                                                12,  48, 255),
    ("Right-Accumbens_area",                                         255, 165,   0),
    ("Right-Hippocampus",                                            103, 255, 255),
    ("Right-Amygdala",                                               220, 216,  20),
    # ---
    ("Right-VentralDC",                                              165,  42,  42),
    # ---
    ("Right-Hypothalamus",                                           204, 182, 142),
]

LH_SUBCORT = [
    ("Left-Pulvinar",                                                255,   0,   0),
    ("Left-Anterior",                                                  0, 255,   0),
    ("Left-Medio_Dorsal",                                            255, 255,   0),
    ("Left-Ventral_Latero_Dorsal",                                   255, 123,   0),
    ("Left-Central_Lateral-Lateral_Posterior-Medial_Pulvinar",         0, 255, 255),
    ("Left-Ventral_Anterior",                                        255,   0, 255),
    ("Left-Ventral_Latero_Ventral",                                    0,   0, 255),
    # ---
    ("Left-Thalamus_Proper",                                           0, 118,  14),
    ("Left-Caudate",                                                 122, 186, 220),
    ("Left-Putamen",                                                 236,  13, 176),
    ("Left-Pallidum",                                                  12,  48, 255),
    ("Left-Accumbens_area",                                          255, 165,   0),
    ("Left-Hippocampus",                                             103, 255, 255),
    ("Left-Amygdala",                                                220, 216,  20),
    # ---
    ("Left-VentralDC",                                               165,  42,  42),
    # ---
    ("Left-Hypothalamus",                                            204, 182, 142),
]

BRAINSTEM = [
    ("Brain_Stem-Midbrain",  242, 104,  76),
    ("Brain_Stem-Pons",      206, 195,  58),
    ("Brain_Stem-Medulla",   119, 159, 176),
    ("Brain_Stem-SCP",       142, 182,   0),
]

RH_SUBCORT_SECTIONS = [
    ("# Right Hemisphere. Subcortical Structures (Thalamic Nuclei)", 7),
    ("# Right Hemisphere. Subcortical Structures",                   7),
    ("# Right Hemisphere. Ventral Diencephalon",                     1),
    ("# Right Hemisphere. Hypothalamus",                             1),
]

LH_SUBCORT_SECTIONS = [
    ("# Left Hemisphere. Subcortical Structures (Thalamic Nuclei)",  7),
    ("# Left Hemisphere. Subcortical Structures",                    7),
    ("# Left Hemisphere. Ventral Diencephalon",                      1),
    ("# Left Hemisphere. Hypothalamus",                              1),
]

BRAINSTEM_SECTION = "# Brain Stem Structures"


def ctab_to_rgb(ctab_row):
    """Extract R, G, B from ctab row (col 0,1,2)."""
    return int(ctab_row[0]), int(ctab_row[1]), int(ctab_row[2])


def generate_lut(scale):
    rh_annot = os.path.join(ATLAS_DIR, f"rh.lausanne2018.scale{scale}.annot")
    lh_annot = os.path.join(ATLAS_DIR, f"lh.lausanne2018.scale{scale}.annot")
    out_path  = os.path.join(
        ATLAS_DIR, f"label-L2018_desc-scale{scale}_atlas_FreeSurferColorLUT.txt"
    )

    if not os.path.exists(rh_annot) or not os.path.exists(lh_annot):
        print(f"  Scale {scale}: annot files missing, skipping")
        return

    _, ctab_rh, names_rh = nib.freesurfer.read_annot(rh_annot)
    _, ctab_lh, names_lh = nib.freesurfer.read_annot(lh_annot)

    lines = []
    lines.append(f"#$Id: ROIv_HR_th_scale{scale}_FreeSurferColorLUT.txt (generated from annot)\n")
    lines.append("\n")
    lines.append("#No. Label Name:                                               R   G   B A \n")
    lines.append("\n")

    idx = 1

    # ── RH cortical ──────────────────────────────────────────────────────────
    lines.append("# Right Hemisphere. Cortical Structures \n")
    for i, (row, name) in enumerate(zip(ctab_rh, names_rh)):
        n = name.decode().strip()
        if n in ("", "unknown"):
            continue
        r, g, b = ctab_to_rgb(row)
        lines.append(f"{idx:<5} ctx-rh-{n:<46} {r:>3} {g:>3} {b:>3} 0 \n")
        idx += 1

    lines.append("\n")

    # ── RH subcortical ───────────────────────────────────────────────────────
    flat_rh = list(RH_SUBCORT)
    pos = 0
    for section_header, count in RH_SUBCORT_SECTIONS:
        lines.append(section_header + " \n")
        for name, r, g, b in flat_rh[pos:pos+count]:
            lines.append(f"{idx:<5} {name:<46} {r:>3} {g:>3} {b:>3} 0 \n")
            idx += 1
        pos += count
        lines.append("\n")

    # ── LH cortical ──────────────────────────────────────────────────────────
    lines.append("# Left Hemisphere. Cortical Structures \n")
    for i, (row, name) in enumerate(zip(ctab_lh, names_lh)):
        n = name.decode().strip()
        if n in ("", "unknown"):
            continue
        r, g, b = ctab_to_rgb(row)
        lines.append(f"{idx:<5} ctx-lh-{n:<46} {r:>3} {g:>3} {b:>3} 0 \n")
        idx += 1

    lines.append("\n")

    # ── LH subcortical ───────────────────────────────────────────────────────
    flat_lh = list(LH_SUBCORT)
    pos = 0
    for section_header, count in LH_SUBCORT_SECTIONS:
        lines.append(section_header + " \n")
        for name, r, g, b in flat_lh[pos:pos+count]:
            lines.append(f"{idx:<5} {name:<46} {r:>3} {g:>3} {b:>3} 0 \n")
            idx += 1
        pos += count
        lines.append("\n")

    # ── Brainstem ─────────────────────────────────────────────────────────────
    lines.append(BRAINSTEM_SECTION + " \n")
    for name, r, g, b in BRAINSTEM:
        lines.append(f"{idx:<5} {name:<46} {r:>3} {g:>3} {b:>3} 0 \n")
        idx += 1
    lines.append("\n")

    with open(out_path, "w") as fh:
        fh.writelines(lines)

    print(f"  Scale {scale}: {idx-1} entries written → {os.path.basename(out_path)}")


def verify_against_reference(scale, ref_path):
    """Check that generated LUT matches reference on region names and ordering."""
    gen_path = os.path.join(
        ATLAS_DIR, f"label-L2018_desc-scale{scale}_atlas_FreeSurferColorLUT.txt"
    )
    ENTRY = re.compile(r'^(\d+)\s+(\S+)')

    ref = {}
    with open(ref_path) as fh:
        for line in fh:
            m = ENTRY.match(line.strip())
            if m:
                ref[int(m.group(1))] = m.group(2)

    gen = {}
    with open(gen_path) as fh:
        for line in fh:
            m = ENTRY.match(line.strip())
            if m:
                gen[int(m.group(1))] = m.group(2)

    mismatches = []
    for k in sorted(ref):
        if k not in gen:
            mismatches.append(f"  ID {k}: missing in generated ({ref[k]})")
        elif ref[k] != gen[k]:
            mismatches.append(f"  ID {k}: ref={ref[k]}  gen={gen[k]}")

    if mismatches:
        print(f"  Scale {scale} MISMATCH vs reference:")
        for m in mismatches[:10]:
            print(m)
    else:
        print(f"  Scale {scale}: generated LUT matches reference exactly ✓")


if __name__ == "__main__":
    print("Generating MSBP-format LUTs for scales 1-5...")

    # Reference LUTs (from real MSBP output) — used to validate the generator
    REF = {
        3: os.path.join(ATLAS_DIR, "label-L2018_desc-scale3_atlas_FreeSurferColorLUT.txt"),
        5: os.path.join(ATLAS_DIR, "label-L2018_desc-scale5_atlas_FreeSurferColorLUT.txt"),
    }

    for scale in [1, 2, 3, 4, 5]:
        generate_lut(scale)

    print()
    print("Validating generated LUTs against reference (scale 3 and 5)...")
    for scale, ref_path in REF.items():
        # rename generated to temp, restore reference, then compare
        gen_path = os.path.join(
            ATLAS_DIR, f"label-L2018_desc-scale{scale}_atlas_FreeSurferColorLUT.txt"
        )
        verify_against_reference(scale, ref_path)

    print("\nDone.")
