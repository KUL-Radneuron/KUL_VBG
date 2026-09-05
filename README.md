# KUL_VBG

**KUL_VBG** ("KULeuven - Virtual Brain Grafting") enables whole-brain parcellation in the presence of large lesions.

Whole-brain parcellation means labeling the brain into parts, gyri, etc. Many parcellation software packages fail in the presence of large brain lesions.

The approach taken here is to:

1. **Extract** the gross brain lesion using a mask
2. **Replace** the brain lesion with normal-looking tissue (hence "virtual brain grafting")
3. **Run** parcellation software like FreeSurfer and/or FastSurfer (which only work well on non-lesioned brains)
4. **Reinsert** the lesion mask into the parcellation

---

## Introduction

For reference, see the paper published in *Neuroimage*, 2021: <https://doi.org/10.1016/j.neuroimage.2021.117731>

> **Which branch to use:** use the `master` branch — it is the updated and more stable version of KUL_VBG. For the version corresponding to the published article, see the branch `Orig_doi.org/10.1016/j.neuroimage.2021.117731`.

---

## Posing the problem and solution

**The problem** — FreeSurfer will not parcellate these brains:

![VBG fig1](figs4readme/fig1.jpg)

**The solution** — the graphical overview of the VBG workflow:

![VBG fig2](figs4readme/fig2.jpg)

---

## Using VBG

VBG has been tested on macOS, WSL2 on Windows 11, Ubuntu 18.0/20.04, Mint 20/21, and CentOS.

VBG can be installed locally after installing the dependencies listed below, or run via a container image — see [**Using VBG via Docker or Apptainer**](#using-vbg-via-docker-or-apptainer) below. Apptainer in particular needs no root at all, which makes it the practical option on a shared/HPC server.

### Dependencies for installing VBG locally (v2.0)

| Dependency | Version | Notes |
|---|---|---|
| **FreeSurfer** | `>= 7.3` (validated on 8.2.0) | **Required, not optional.** VBG 2.0 uses `mri_synthstrip` for brain extraction (`-B 1`), `mri_synthseg` (`-P 4`), and `segment_subregions` / `mri_segment_hypothalamic_subunits` for the multi-atlas block (`-M`). None of these exist in FreeSurfer 6. |
| **ANTs** | v2.4.4 | plus the ANTsX scripts |
| **FSL** | v6.0.7 | only `fslmaths`, `fslstats`, `fslreorient2std`, `fslswapdim`, `fslorient` and `imcp` are actually used |
| **MRtrix3** | v3.0.4+ | command-line tools only; `mrview` is not needed |
| **FastSurfer** | — | only required for `-P 2` / `-P 3` |
| **Python 3** | — | with `numpy`, `scipy`, `nibabel`, `matplotlib` (for `KUL_VBG_QC.py` and `KUL_lesion_overlap.py`) |

> **Note:** HD-BET is no longer a dependency. VBG 2.0 replaced both HD-BET code paths with SynthStrip (`-B 1`, default) and ANTs-BET (`-B 2`).

For help setting up your environment with the various neuroimaging packages, see [KUL_Linux_Installation](https://github.com/treanus/KUL_Linux_Installation.git).

---

## Examples

### Inputs

**Obligatory:**

| Flag | Meaning |
|---|---|
| `-S` | Subject/participant name in BIDS convention, without the leading `sub-` |
| `-a` | A NIfTI-format T1WI of the subject |
| `-l` | Binary lesion mask (lesion = 1, background = 0), integer NIfTI format |
| `-z` | Lesion mask space. **N.B.** the specified lesion mask must have the same dimensions and transform as the input T1WI. |

**Optional:**

| Flag | Meaning |
|---|---|
| `-m`, `-o` | Location of intermediate processing and output folders |
| `-n` | Number of parallel workers used |
| `-t` | Type of filling (default = uVBG; add `-t` to activate bVBG) |
| `-p` | Age group of participant (default = adult; add `-p` for pediatric-friendly mode) |
| `-P` | Run parcellation *after* lesion filling finishes: `1` = FreeSurfer, `2` = FastSurfer, `3` = FastSurfer/FreeSurfer hybrid |
| `-v` | Verbose mode |

### Example commands

Unilateral VBG approach, HD-BET-equivalent brain extraction, BIDS input with a single session, FreeSurfer for parcellation:

```bash
KUL_VBG.sh -p pat001 -b -n 6 -l /fullpath/lesion_T1w.nii.gz -z T1 -o /fullpath/output -B 1 -P 1 -v
```

Bilateral VBG approach, non-BIDS input, FastSurfer for parcellation:

```bash
KUL_VBG.sh -p pat001 -a /fullpath/sub-PT_T1w.nii.gz -n 6 -l /fullpath/lesion_T1w.nii.gz -z T1 -o /fullpath/output -t -B 1 -P 2 -v
```

### Purpose

The purpose of this workflow is to generate a lesion-filled image, with healthy-looking synthetic tissue in place of the lesion — essentially excising the lesion and grafting over the resulting defect in the T1 MR image space.

### Required arguments

| Flag | Meaning |
|---|---|
| `-S` | BIDS participant name (anonymised name of the subject without the `sub-` prefix) |
| `-b` | Data is in BIDS |
| `-l` | Full path and file name to the lesion mask file per session |
| `-z` | Space of the lesion mask used (only T1 supported in this version) |
| `-a` | Input precontrast T1WIs |

### Optional arguments

| Flag | Meaning |
|---|---|
| `-s` | Session (of the participant) |
| `-t` | Use the VBG template to derive the fill patch (if used, template tissue is used alongside native tissue to create the donor brain) |
| `-E` | Treat as an extra-axial lesion (skip VBG bulk, fill lesion patch with 0s, run FS and subsequent steps) |
| `-B` | Brain extraction method: `1` = HD-BET, `2` = ANTs-BET (default if not set) |
| `-P` | Run parcellation: `1` = FreeSurfer recon-all, `2` = FastSurfer |
| `-p` | Pediatric patients — use pediatric template (`NKI_under_10` in MNI) |
| `-m` | Full path to intermediate output dir (default: `./VBG_out/proc_VBG`) |
| `-o` | Full path to output dir (default: `./VBG_out/output_VBG`) |
| `-n` | Number of CPUs for parallelisation (default: 6) |
| `-v` | Show output from MRtrix commands |
| `-h` | Prints help menu |

### Notes

- Flags `-b` and `-a` are **mutually exclusive** — if your data is in BIDS use `-b`; if not, specify the exact path and name for the patient's `T1.nii.gz`.
- You need a high-resolution T1WI and a lesion mask in the same space for VBG to run.
- If you end up with an empty image, you likely have a mismatch between the T1 and lesion mask.
- The lesion mask can be generated with any lesion segmentation tool.
- The lesion mask needs to be specific to the lesion, with voxel value `1` encoding the lesion and `0` for healthy tissue.

### Installation instructions

1. Clone this repository and add the installation directory to your `PATH` in the Bash shell.
2. Ensure all dependencies are met. FastSurfer is only required if you will use it for parcellation (i.e. with `-P 2` or `-P 3`).

---

## Using VBG via Docker or Apptainer

Container images bundle FreeSurfer 8.2.0, FSL, ANTs, MRtrix3, FastSurfer and VBG itself, so there is nothing to install beyond the container runtime.

Full documentation, build instructions and troubleshooting are in [`Docker/README.md`](Docker/README.md).

> **You must supply your own FreeSurfer license.** It is deliberately not included in the image — the license is issued per user and cannot be redistributed. Register for free at <https://surfer.nmr.mgh.harvard.edu/registration.html>, then bind it in at `/licence/license.txt` as shown below.

### Apptainer (no root required — use this on a shared or HPC server)

```bash
apptainer run --nv \
    -B /path/to/license.txt:/licence/license.txt \
    -B "$PWD":/data \
    KUL_VBG_2.0.sif \
    KUL_VBG.sh -S PAT001 -a /data/sub-PAT001_T1w.nii.gz \
               -l /data/lesion.nii.gz -z T1 \
               -o /data/VBG_out -B 1 -P 1 -n 8 -v
```

### Docker

```bash
docker run --gpus all --rm -it \
    -v /path/to/license.txt:/licence/license.txt:ro \
    -v "$PWD":/data \
    kul_vbg:2.0 \
    KUL_VBG.sh -S PAT001 -a /data/sub-PAT001_T1w.nii.gz \
               -l /data/lesion.nii.gz -z T1 \
               -o /data/VBG_out -B 1 -P 1 -n 8 -v
```

Drop `--nv` / `--gpus all` if there is no GPU — everything still works, and only FastSurfer (`-P 2`/`-P 3`) is slower.

> **Note:** all paths passed to `KUL_VBG.sh` must be paths **inside** the container (`/data/...`), not host paths.

### Building the images

```bash
cd Docker
./build.sh --sif     # builds the Docker image, then converts it to a .sif
```

The resulting `.sif` is a single self-contained file: copy it to the server and run it there — no root and no installation required.
