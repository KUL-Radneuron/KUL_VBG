#!/usr/bin/env python3
"""
make_readable_lut.py

Convert a Lausanne 2018 FreeSurfer LUT to:
  - short, human-readable label names (R.SupFront_3 instead of ctx-rh-superiorfrontal_3)
  - lobar-colour-coded cortical parcels (consistent hue family per lobe)
  - subcortical entries kept with their original colours

Usage:
  python3 make_readable_lut.py input_lut.txt output_lut.txt
"""

import sys
import re
import colorsys

# ── Short display names ───────────────────────────────────────────────────────
ABBREV = {
    "lateralorbitofrontal":     "LatOrbFront",
    "parsorbitalis":            "ParsOrbit",
    "frontalpole":              "FrontPole",
    "medialorbitofrontal":      "MedOrbFront",
    "parstriangularis":         "ParsTri",
    "parsopercularis":          "ParsOper",
    "rostralmiddlefrontal":     "RMidFront",
    "superiorfrontal":          "SupFront",
    "caudalmiddlefrontal":      "CMidFront",
    "precentral":               "Precentral",
    "paracentral":              "Paracentral",
    "rostralanteriorcingulate": "rACC",
    "caudalanteriorcingulate":  "cACC",
    "posteriorcingulate":       "PostCing",
    "isthmuscingulate":         "IsthmusCing",
    "postcentral":              "Postcentral",
    "supramarginal":            "Supramarg",
    "superiorparietal":         "SupPar",
    "inferiorparietal":         "InfPar",
    "precuneus":                "Precuneus",
    "cuneus":                   "Cuneus",
    "pericalcarine":            "Pericalc",
    "lateraloccipital":         "LatOccip",
    "lingual":                  "Lingual",
    "fusiform":                 "Fusiform",
    "parahippocampal":          "ParaHipp",
    "entorhinal":               "Entorh",
    "temporalpole":             "TempPole",
    "inferiortemporal":         "InfTemp",
    "middletemporal":           "MidTemp",
    "bankssts":                 "BankSTS",
    "superiortemporal":         "SupTemp",
    "transversetemporal":       "TransTemp",
    "insula":                   "Insula",
}

# ── Lobar colour families: (hue_degrees, sat_max, val_max) ───────────────────
# Hue families chosen to be maximally distinct from each other:
#   Frontal=blue  Motor=red  Cingulate=purple  Parietal=green
#   Occipital=yellow-green  Temporal=orange  Insula=teal
LOBE_PARAMS = {
    "frontal":    (225, 0.90, 0.92),
    "motor":      (  5, 0.90, 0.88),
    "cingulate":  (295, 0.80, 0.85),
    "parietal":   (135, 0.85, 0.82),
    "occipital":  ( 72, 0.88, 0.88),
    "temporal":   ( 28, 0.88, 0.88),
    "insula":     (180, 0.82, 0.82),
}

REGION_LOBE = {
    "lateralorbitofrontal":     "frontal",
    "parsorbitalis":            "frontal",
    "frontalpole":              "frontal",
    "medialorbitofrontal":      "frontal",
    "parstriangularis":         "frontal",
    "parsopercularis":          "frontal",
    "rostralmiddlefrontal":     "frontal",
    "superiorfrontal":          "frontal",
    "caudalmiddlefrontal":      "frontal",
    "precentral":               "motor",
    "paracentral":              "motor",
    "postcentral":              "motor",
    "rostralanteriorcingulate": "cingulate",
    "caudalanteriorcingulate":  "cingulate",
    "posteriorcingulate":       "cingulate",
    "isthmuscingulate":         "cingulate",
    "supramarginal":            "parietal",
    "superiorparietal":         "parietal",
    "inferiorparietal":         "parietal",
    "precuneus":                "parietal",
    "cuneus":                   "occipital",
    "pericalcarine":            "occipital",
    "lateraloccipital":         "occipital",
    "lingual":                  "occipital",
    "fusiform":                 "temporal",
    "parahippocampal":          "temporal",
    "entorhinal":               "temporal",
    "temporalpole":             "temporal",
    "inferiortemporal":         "temporal",
    "middletemporal":           "temporal",
    "bankssts":                 "temporal",
    "superiortemporal":         "temporal",
    "transversetemporal":       "temporal",
    "insula":                   "insula",
}

CTX_RE = re.compile(r'^ctx-(lh|rh)-([a-zA-Z]+)_(\d+)$')


def parcel_rgb(region, parcel_num, hemi):
    """Golden-ratio hue/sat/val spread within a lobe colour band.
    Region name is hashed into the seed so that ParsOrbit_1 and SupFront_1
    land at different points even though both have parcel_num=1.
    """
    lobe = REGION_LOBE.get(region, "frontal")
    hue_deg, sat_max, val_max = LOBE_PARAMS[lobe]

    # Region-specific offset (0-99) derived from region name
    region_seed = sum(ord(c) * (i + 1) for i, c in enumerate(region)) % 100
    combined = region_seed * 50 + parcel_num  # unique per (region, parcel_num)

    golden = 0.6180339887498949
    t  = (combined * golden) % 1.0
    t2 = (combined * 0.38196601125) % 1.0

    hue = ((hue_deg + (t - 0.5) * 36) % 360) / 360.0  # ±18 deg jitter
    sat = 0.55 + t  * (sat_max - 0.55)
    val = 0.62 + t2 * (val_max - 0.62)

    # right hemisphere slightly brighter for hemisphere disambiguation
    val = min(1.0, val + 0.04) if hemi == "rh" else max(0.50, val - 0.04)

    r, g, b = colorsys.hsv_to_rgb(hue, sat, val)
    return int(r * 255 + 0.5), int(g * 255 + 0.5), int(b * 255 + 0.5)


def short_name(label):
    m = CTX_RE.match(label)
    if not m:
        return label
    hemi, region, num = m.group(1), m.group(2), m.group(3)
    prefix = "R" if hemi == "rh" else "L"
    abbr = ABBREV.get(region, region[:12])
    return f"{prefix}.{abbr}_{num}"


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input_lut.txt output_lut.txt")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]
    ENTRY_RE = re.compile(r'^(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)')

    out_lines = [
        f"#$Id: {out_path} — human-readable lobar-colour LUT\n",
        "\n",
        "#No. Label Name                                      R   G   B  A\n",
        "\n",
    ]

    with open(in_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            if line.startswith("#"):
                # drop original header line, keep section comments
                if not line.startswith("#$Id"):
                    out_lines.append(line + "\n")
                continue

            if not line.strip():
                out_lines.append("\n")
                continue

            m = ENTRY_RE.match(line.strip())
            if not m:
                out_lines.append(line + "\n")
                continue

            idx = int(m.group(1))
            label = m.group(2)
            orig_r, orig_g, orig_b, orig_a = (
                int(m.group(3)), int(m.group(4)),
                int(m.group(5)), int(m.group(6)),
            )

            cx = CTX_RE.match(label)
            if cx:
                hemi, region, num = cx.group(1), cx.group(2), int(cx.group(3))
                name = short_name(label)
                r, g, b = parcel_rgb(region, num, hemi)
            else:
                name = label
                r, g, b = orig_r, orig_g, orig_b

            out_lines.append(
                f"{idx:<5} {name:<46} {r:>3} {g:>3} {b:>3}  {orig_a}\n"
            )

    with open(out_path, "w") as fh:
        fh.writelines(out_lines)

    print(f"Written {out_path}")


if __name__ == "__main__":
    main()
