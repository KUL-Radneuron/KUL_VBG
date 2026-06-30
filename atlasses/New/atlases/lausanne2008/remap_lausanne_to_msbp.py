#!/usr/bin/env python3
"""
remap_lausanne_to_msbp.py

After mri_aparc2aseg writes a Lausanne parcellation with FS-native label IDs
(1000 + ctab_index for RH, 2000 + ctab_index for LH), this script remaps
every voxel to the MSBP sequential ID scheme so downstream tractography
scripts that hardcode parcel IDs continue to work.

Usage:
  python3 remap_lausanne_to_msbp.py \\
      --input   <raw_aparc2aseg.mgz>      \\
      --lh_annot <lh.lausanne2018.scaleN.annot> \\
      --rh_annot <rh.lausanne2018.scaleN.annot> \\
      --lut     <label-L2018_desc-scaleN_atlas_FreeSurferColorLUT.txt> \\
      --output  <remapped_parcellation.mgz>
"""

import sys
import argparse
import re
import numpy as np
import nibabel as nib


def load_annot_ctab(annot_path):
    """Return {ctab_index: short_name} (index 0 = 'unknown', skipped)."""
    _, _, names = nib.freesurfer.read_annot(annot_path)
    return {i: n.decode().strip() for i, n in enumerate(names)}


def load_msbp_lut(lut_path):
    """Return {full_label_name: msbp_id}.
    full_label_name examples: 'ctx-rh-precentral_4', 'Left-Hippocampus'
    Also builds a normalised lowercase lookup for subcortical matching.
    """
    id_map = {}
    norm_map = {}   # normalised name → msbp_id (for subcortical fuzzy matching)
    ENTRY = re.compile(r'^(\d+)\s+(\S+)')
    with open(lut_path) as fh:
        for line in fh:
            m = ENTRY.match(line.strip())
            if m:
                msbp_id = int(m.group(1))
                name = m.group(2)
                id_map[name] = msbp_id
                norm_map[name.lower().replace('-', '_')] = msbp_id
    return id_map, norm_map


def build_volume_to_msbp(ctab_lh, ctab_rh, msbp_map, norm_map):
    """Build {fs_volume_label: msbp_id} for all cortical entries.

    FS convention:
      RH cortical voxel value = 1000 + ctab_index
      LH cortical voxel value = 2000 + ctab_index
    """
    vol2msbp = {}

    # mri_aparc2aseg with custom annots uses 1000+idx for LH, 2000+idx for RH
    # (opposite to the standard DK atlas convention — confirmed empirically)
    for idx, short_name in ctab_lh.items():
        if idx == 0 or not short_name:
            continue
        fs_id = 1000 + idx
        full_name = f"ctx-lh-{short_name}"
        if full_name in msbp_map:
            vol2msbp[fs_id] = msbp_map[full_name]

    for idx, short_name in ctab_rh.items():
        if idx == 0 or not short_name:
            continue
        fs_id = 2000 + idx
        full_name = f"ctx-rh-{short_name}"
        if full_name in msbp_map:
            vol2msbp[fs_id] = msbp_map[full_name]

    return vol2msbp


# Standard FS aseg label IDs → base names (hyphens, as in FS)
# Used to match subcortical entries in the MSBP LUT
FS_ASEG = {
    16:  "Brain-Stem",
    17:  "Left-Hippocampus",
    18:  "Left-Amygdala",
    11:  "Left-Caudate",
    12:  "Left-Putamen",
    13:  "Left-Pallidum",
    10:  "Left-Thalamus-Proper",
    26:  "Left-Accumbens-area",
    28:  "Left-VentralDC",
    15:  "Left-Hypothalamus",      # not standard aseg but MSBP includes it
    53:  "Right-Hippocampus",
    54:  "Right-Amygdala",
    50:  "Right-Caudate",
    51:  "Right-Putamen",
    52:  "Right-Pallidum",
    49:  "Right-Thalamus-Proper",
    58:  "Right-Accumbens-area",
    60:  "Right-VentralDC",
}

# MSBP uses underscores in some names where FS uses hyphens;
# build a normalised lookup that strips such differences
def add_subcortical(vol2msbp, norm_map):
    for fs_id, fs_name in FS_ASEG.items():
        norm = fs_name.lower().replace('-', '_')
        # direct match
        if norm in norm_map:
            vol2msbp[fs_id] = norm_map[norm]
            continue
        # try dropping _proper / _area suffix variants
        for key in norm_map:
            if key.startswith(norm[:12]):
                vol2msbp[fs_id] = norm_map[key]
                break


def remap_volume(vol_data, vol2msbp):
    """Apply vol2msbp lookup; unmapped labels become 0."""
    out = np.zeros_like(vol_data, dtype=np.int32)
    unique = np.unique(vol_data)
    for uid in unique:
        if uid == 0:
            continue
        target = vol2msbp.get(int(uid), 0)
        if target:
            out[vol_data == uid] = target
        # else: label not in LUT → stays 0 (background / not in atlas)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',    required=True, help='raw mri_aparc2aseg output')
    ap.add_argument('--lh_annot', required=True, help='lh.lausanne2018.scaleN.annot')
    ap.add_argument('--rh_annot', required=True, help='rh.lausanne2018.scaleN.annot')
    ap.add_argument('--lut',      required=True, help='MSBP FreeSurferColorLUT.txt')
    ap.add_argument('--output',   required=True, help='remapped parcellation .mgz')
    args = ap.parse_args()

    print("Loading annot ctabs...")
    ctab_lh = load_annot_ctab(args.lh_annot)
    ctab_rh = load_annot_ctab(args.rh_annot)

    print("Loading MSBP LUT...")
    msbp_map, norm_map = load_msbp_lut(args.lut)

    print("Building volume-ID → MSBP-ID mapping...")
    vol2msbp = build_volume_to_msbp(ctab_lh, ctab_rh, msbp_map, norm_map)
    add_subcortical(vol2msbp, norm_map)

    print("Loading raw parcellation volume...")
    img = nib.load(args.input)
    vol = np.asarray(img.dataobj, dtype=np.int32)

    # Diagnostic: show what unique labels are present and whether they map
    unique_labels = np.unique(vol)
    n_mapped   = sum(1 for u in unique_labels if u != 0 and int(u) in vol2msbp)
    n_unmapped = sum(1 for u in unique_labels if u != 0 and int(u) not in vol2msbp)
    print(f"  Unique non-zero labels: {len(unique_labels)-1}")
    print(f"  Mapped to MSBP IDs:     {n_mapped}")
    print(f"  Unmapped (→ 0):         {n_unmapped}")
    if n_unmapped > 0:
        unmapped = [int(u) for u in unique_labels if u != 0 and int(u) not in vol2msbp]
        print(f"  Unmapped label values:  {unmapped[:20]}")

    print("Remapping...")
    out_vol = remap_volume(vol, vol2msbp)

    out_img = nib.MGHImage(out_vol, img.affine, img.header)
    nib.save(out_img, args.output)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
