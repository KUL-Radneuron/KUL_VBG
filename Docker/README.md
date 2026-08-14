# KUL_VBG 2.0 — Docker and Apptainer

Container images for [KUL_VBG](https://github.com/KUL-Radneuron/KUL_VBG), the
Virtual Brain Grafting workflow: whole-brain parcellation in the presence of
large lesions.

If you use this in research, please cite:

> Radwan AM, Emsell L, Blommaert J, Zhylka A, Kovacs S, Theys T, Sollmann N,
> Dupont P, Sunaert S. **Virtual brain grafting: Enabling whole brain
> parcellation in the presence of large lesions.** *NeuroImage* 2021;229:117731.
> <https://doi.org/10.1016/j.neuroimage.2021.117731>

---

## Quick start

### Run with Apptainer (rootless servers, HPC)

```bash
apptainer run --nv \
    -B /path/to/license.txt:/licence/license.txt \
    -B "$PWD":/data \
    KUL_VBG_2.0.sif \
    KUL_VBG.sh -S PAT001 -a /data/sub-PAT001_T1w.nii.gz \
               -l /data/lesion.nii.gz -z T1 \
               -o /data/VBG_out -B 1 -P 1 -n 8 -v
```

### Run with Docker

```bash
docker run --gpus all --rm -it \
    -v /path/to/license.txt:/licence/license.txt:ro \
    -v "$PWD":/data \
    kul_vbg:2.0 \
    KUL_VBG.sh -S PAT001 -a /data/sub-PAT001_T1w.nii.gz \
               -l /data/lesion.nii.gz -z T1 \
               -o /data/VBG_out -B 1 -P 1 -n 8 -v
```

Drop `--nv` / `--gpus all` if there is no GPU. Everything still works;
FastSurfer (`-P 2`/`-P 3`) just runs on CPU.

> **Paths are container paths.** Everything you pass to `KUL_VBG.sh` must be the
> path *inside* the container (`/data/...`), not the host path.

---

## You must supply a FreeSurfer licence

The licence is **not** included in the image, and never will be: FreeSurfer's
licence is issued per user and cannot be redistributed inside a shared image.

1. Register (free): <https://surfer.nmr.mgh.harvard.edu/registration.html>
2. Bind the resulting `license.txt` at **`/licence/license.txt`**:
   - Docker: `-v /path/to/license.txt:/licence/license.txt:ro`
   - Apptainer: `-B /path/to/license.txt:/licence/license.txt`

Inside the image `$FREESURFER_HOME/license.txt` is a symlink pointing at that
path, which is what makes this work on a read-only Apptainer SIF — nothing is
ever copied or written inside the image. If you forget it, the entrypoint stops
immediately with instructions rather than failing deep inside `recon-all`.

---

## Building

Everything is driven by `build.sh` in this directory.

```bash
cd Docker

./build.sh              # Docker image only
./build.sh --sif        # Docker image, then convert to KUL_VBG_2.0.sif
./build.sh --help       # all options
```

The Docker image is the single source of truth; the SIF is produced from it, so
the two can never drift apart.

**Expect 1–3 hours and ~40 GB of scratch on a first build** — ANTs and MRtrix3
are compiled from source, and FreeSurfer 8.2.0 plus FastSurfer's checkpoints are
downloaded. Later builds reuse the layer cache.

`build.sh` smoke-tests the image afterwards (every expected binary on `PATH`,
Python QC imports, atlas data present), so a broken build fails at build time
rather than an hour into a patient run.

> **Install `docker buildx` first if you can.** The ANTs, MRtrix3 and FSL stages
> are independent, so BuildKit builds them **in parallel** and skips stages a
> given target does not need. The deprecated legacy builder does neither — it
> walks every stage sequentially, so even `--target fsl-builder` compiles ANTs
> first. The build succeeds either way, just markedly slower.
>
> ```bash
> sudo apt-get install -y docker-buildx
> ```
>
> `build.sh` detects which builder is available and tells you which it is using.

### Choosing a CUDA backend

This is the one build option worth thinking about. FastSurfer's PyTorch backend
**cannot be detected automatically**, because `docker build` has no GPU.

| `--torch-backend` | Use when |
|---|---|
| `cu126` *(default)* | Almost everything. Driver ≥ 525, via CUDA minor-version compatibility. |
| `cu128` | Blackwell-era GPUs on a recent driver. |
| `cu118` | Very old drivers. |
| `cpu` | No GPU anywhere — much smaller image. |

Run `nvidia-smi` on the **target server**: the "CUDA Version" it reports is the
highest your driver supports. Pick a backend at or below that.

### FSL: subset vs full

By default the image installs an FSL **subset** — 721 MB instead of ~11 GB.

This is measured, not guessed. The complete set of FSL calls in the KUL_VBG
source tree is:

| Binary | Call sites | Conda package |
|---|---|---|
| `fslmaths` | 114 | `fsl-avwutils` |
| `fslstats` | 15 | `fsl-avwutils` |
| `fslreorient2std` | 6 | `fsl-avwutils` |
| `fslswapdim`, `fslorient` | 4 | `fsl-avwutils` |
| `imcp` | 1 | `fslpy` |
| `MNI152_T1_1mm_brain.nii.gz` (QC) | — | `fsl-data_standard` |

VBG does brain extraction with SynthStrip/ANTs and tissue segmentation with ANTs
Atropos, so it never calls FSL's own `bet`, `fast` or `first`.

Crucially this is a **conda package** subset, not hand-copied binaries — conda
resolves the shared-library closure automatically, which is where subset
installs normally break. If some rarely-hit path ever calls an FSL tool outside
this set, it fails loudly with "command not found" at that call, never silently
or with wrong numbers.

To install the full distribution instead:

```bash
./build.sh --fsl full
```

---

## Deploying to a server where you have no root

This is the intended workflow, and it needs no privileges on the server at all.

```bash
# 1. On your own machine (Docker + Apptainer available):
cd Docker
./build.sh --sif

# 2. Copy the single self-contained file across:
scp KUL_VBG_2.0.sif user@server:~/containers/

# 3. On the server — no root, no install step:
apptainer run --nv \
    -B ~/license.txt:/licence/license.txt \
    -B /data/study:/data \
    ~/containers/KUL_VBG_2.0.sif \
    KUL_VBG.sh -S PAT001 -a /data/T1w.nii.gz -l /data/lesion.nii.gz \
               -z T1 -o /data/VBG_out -B 1 -P 1 -n 8 -v
```

The SIF is one file. There is nothing to install, and no root is required to
run it.

### Notes for shared servers

- **`--cleanenv` is worth adding.** Apptainer merges the host environment into
  the container by default. `entrypoint.sh` already pins the variables that
  matter, but on a neuroimaging server where users export their own `FSLDIR` /
  `FREESURFER_HOME`, `--cleanenv` removes a whole class of surprises.
- **`/tmp` must be writable.** KUL_VBG works around a confirmed FreeSurfer 8.2.0
  `mris_register` buffer overflow by routing `SUBJECTS_DIR` through a short
  symlink under `/tmp`. If your site mounts `/tmp` read-only or with a tight
  quota, bind over it: `-B /scratch/$USER/tmp:/tmp`. The entrypoint checks this
  and fails fast with that suggestion.
- **`--nv` is required for GPU.** Without it the host driver is not injected and
  PyTorch sees no GPU *even on a GPU node* — a silent, very slow failure mode.
- **Set `APPTAINER_TMPDIR`** somewhere with room if `/tmp` is small; the
  conversion step needs roughly twice the image size.

### Testing local KUL_VBG changes without rebuilding

The image clones KUL_VBG at a pinned commit, so it always contains committed
code. While developing, bind your working tree over `/opt/KUL_VBG` instead of
rebuilding — the scripts are read from there at run time:

```bash
apptainer run \
    -B /path/to/your/KUL_VBG:/opt/KUL_VBG \
    -B ~/license.txt:/licence/license.txt \
    -B "$PWD":/data \
    KUL_VBG_2.0.sif KUL_VBG.sh -h
```

This picks up uncommitted edits instantly and keeps all the dependencies from
the image. Only rebuild when you want to bake a new pinned commit in:

```bash
./build.sh --sif --build-arg VBG_COMMIT=<sha>
```

### Interactive shell

```bash
apptainer shell -B ~/license.txt:/licence/license.txt -B "$PWD":/data KUL_VBG_2.0.sif
```

Useful for `KUL_synth_pats_4VBG.sh`, `KUL_VBG_cook_template.sh` and
`KUL_VBG_multiparc.sh`, which are all present on `PATH`.

---

## What is in the image

| Component | Version | Notes |
|---|---|---|
| FreeSurfer | 8.2.0 | Installed as a real `.deb` so its ITK dependencies resolve |
| FSL | 6.0.7.23 | Subset by default (see above) |
| ANTs | 2.4.4.post20 (`40ee2d22`) | Built with **gcc-12**, not the Ubuntu 24.04 default |
| MRtrix3 | KUL fork @ `5a643594` | **CLI only** — no Qt, no `mrview` |
| FastSurfer | pinned commit | uv `.venv`; checkpoints baked in |
| KUL_VBG | 2.0 (`145ba69`) | Cloned at a pinned commit |
| Python | numpy, scipy, nibabel, matplotlib | For `KUL_VBG_QC.py`, `KUL_lesion_overlap.py` |

Version pins mirror
[`KUL_Linux_setup/setup_environment.sh`](https://github.com/treanus/KUL_Linux_Installation),
i.e. the versions the pipeline is actually validated against on bare metal.
**Keep the two in sync when bumping either.**

Three of these carry non-obvious reasoning worth preserving:

- **ANTs must be built with gcc-12.** ANTs at this commit vendors a mid-2023
  ITKv5 snapshot whose generated enum headers rely on `<cstdint>` arriving
  transitively. GCC 13 (the Ubuntu 24.04 default) stopped doing that, so the ITK
  sub-build dies with `'uint8_t' was not declared`. This is a documented
  ITK/GCC13 incompatibility ([ITK #4607](https://github.com/InsightSoftwareConsortium/ITK/issues/4607)),
  not something specific to this image.
- **FreeSurfer is installed via `dpkg -i` + `apt-get install -f`, not
  `dpkg-deb -x`.** A payload-only extraction skips the `.deb`'s declared runtime
  dependencies and breaks ITK-linked tools such as `freeview`
  (`libITKCommon-5.3.so.1: cannot open shared object file`) — ITK is not bundled
  inside FreeSurfer's tree, only referenced as a dependency apt must supply.
- **FastSurfer uses `uv pip compile | uv pip sync`, not a bare
  `uv pip sync requirements.txt`.** `requirements.txt` deliberately excludes the
  `nvidia-*` CUDA runtime packages and expects the compile step to pull them in
  transitively. `sync` alone installs only what is literally listed, silently
  omitting the runtime libs and leaving torch unable to import
  (`libcudart.so.*: cannot open shared object file`).

FastSurfer's checkpoints are downloaded **at build time**. This is the single
most important difference between an image that works on an HPC node and one
that does not: otherwise FastSurfer fetches them on first run, which fails on a
compute node with no outbound network and tries to write into a read-only SIF.

---

## Upgrading from the v1.x image

The old `KUL_VBG_Dockerfile` will not produce a working VBG 2.0 image. What
changed, and why:

| Old | Now |
|---|---|
| FreeSurfer 6.0.0 | **8.2.0** — VBG 2.0 needs FS ≥ 7.3 for `mri_synthstrip` (`-B 1`), `mri_synthseg` (`-P 4`), `segment_subregions` and `mri_segment_hypothalamic_subunits` (`-M`). None exist in FS 6. |
| HD-BET + model weights | **Removed.** VBG 2.0 dropped both HD-BET code paths in favour of SynthStrip (`-B 1`) / ANTs-BET (`-B 2`). |
| Docker-in-Docker, NVIDIA Container Toolkit | **Removed.** Always questionable inside an image; unusable under Apptainer. |
| CUDA toolkit 12.5 (~4 GB) | **Removed.** PyTorch wheels ship their own CUDA runtime; the driver is injected at run time. Nothing here compiles CUDA code. |
| `COPY ANTs` from the build host | Built from a pinned commit — reproducible from the Dockerfile alone. |
| MRtrix3 `dev_sans_cmake` | That branch no longer exists; upstream deleted `./configure && ./build` in Oct 2023. Now a CMake build. |
| `conda env create -f env/fastsurfer.yml` | That file is an **empty placeholder** in current FastSurfer — the old recipe silently produced a broken install. Now uv + `.venv`. |
| Entrypoint `cp`s the licence into the image | Symlink + bind, so it works on a read-only SIF. |
| Linux Mint 21 base | `ubuntu:24.04`, matching the FreeSurfer 8.2.0 `ubuntu24` build. |

---

## Troubleshooting

**"no FreeSurfer licence found"** — bind your licence at
`/licence/license.txt`. See above. Note the image ships an *empty* placeholder
there, and an empty file is rejected on purpose.

**FastSurfer is extremely slow / `torch.cuda.is_available()` is False** — you
almost certainly omitted `--nv` (Apptainer) or `--gpus all` (Docker). The
entrypoint prints a GPU line at startup; check what it says. If it reports a GPU
and torch still disagrees, your `--torch-backend` is likely newer than the host
driver supports — rebuild with a lower one.

**`/tmp is not writable`** — bind a writable directory over it:
`-B /scratch/$USER/tmp:/tmp`.

**Build fails in the ANTs stage with `'uint8_t' was not declared`** — the gcc-12
pin is not taking effect. Check that `gcc-12`/`g++-12` installed correctly in
the `ants-builder` stage.

**Build fails during FastSurfer's `uv pip sync`** — your `--torch-backend` has
no published wheel for the torch version FastSurfer pins. Check
`https://download.pytorch.org/whl/<backend>` and pick one that does.

**Apptainer conversion fails with `archive/tar: invalid tar header`** — this one
message has two completely different causes, and neither means your image is bad.

*Cause 1 — BuildKit attestations.* By default BuildKit attaches provenance/SBOM
attestations, which forces the result into a manifest *list* containing the real
image plus an attestation manifest whose "layer" is an in-toto JSON blob rather
than a tar. Apptainer walks the list and tries to untar the JSON. `build.sh`
passes `--provenance=false --sbom=false`; if you call `docker build` by hand,
pass them yourself. Affects buildx only.

*Cause 2 — apptainer cannot unpack the large FreeSurfer layer (the usual one).*
apptainer 1.5.3 fails partway through this image's ~11 GB FreeSurfer layer,
dying at tar entry 33,830 of 107,050 on an ordinary 1.1 MB file. The layer is
**not** damaged: its blob hashes to exactly the digest it is stored under, and
GNU tar reads all 107,050 entries with exit 0. It is an apptainer-side defect.

`build.sh` handles this automatically. It first tries the normal OCI route, and
on failure falls back to flattening the image so apptainer never parses a layer
at all:

1. `docker export` → one uncompressed tar, no layers, no gzip
2. GNU tar extracts it (under `unshare --map-auto`, so ownership is preserved)
3. apptainer squashes the resulting directory, with `ENV`/`ENTRYPOINT` restored
   through a generated `.def` file (`docker export` carries only the filesystem)

The fast path is attempted first, so this self-heals if a future apptainer
release fixes the bug. The fallback needs roughly 2× the image size in scratch
next to the output SIF.

*Cause 3 — an intermittently corrupt layer export.* Observed here on
docker 29.1.3: exporting this image's ~11 GB FreeSurfer layer sometimes produced
a blob of the correct length but wrong bytes, so it no longer hashed to the
digest it was stored under. Re-exporting the identical image produced a clean
blob, so the stored layer was never damaged — only the export, and only
sometimes. Smaller layers were unaffected every time.

`build.sh` now defends against this: it exports to an archive, checks every blob
against its own digest (content-addressed blobs make this exact and free — the
name *is* the checksum), and retries a bad export up to 3 times before giving
up. Verification takes ~2 minutes and replaces discovering the problem 20
minutes into a conversion. Cost: the archive needs ~16 GB of scratch alongside
the image.

If it fails all 3 attempts, suspect hardware rather than Docker — intermittent
corruption of only the largest object, with nothing in `dmesg`, is a classic
non-ECC RAM bit-flip signature. Run `memtest86+`.

**`command not found` for an FSL tool** — you have hit a code path outside the
measured subset. Rebuild with `./build.sh --fsl full` and please open an issue
so the subset list can be corrected.
