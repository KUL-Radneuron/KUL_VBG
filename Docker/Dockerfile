# syntax=docker/dockerfile:1.7
###############################################################################
# KUL_VBG 2.0 — container image
#
# Ahmed M. Radwan / Stefan Sunaert, KU Leuven — radwanphd@gmail.com
#
# Built to be run BOTH as a Docker image and, after conversion, as an Apptainer
# SIF on a rootless multi-user server. See README.md in this directory.
#
# Version pins below mirror KUL_Linux_setup/setup_environment.sh — i.e. the
# versions this pipeline is actually validated against on bare metal, not
# necessarily the newest upstream. Keep the two in sync when bumping either.
#
# ---------------------------------------------------------------------------
# WHAT CHANGED vs the old KUL_VBG_Dockerfile (v1.x), and why
# ---------------------------------------------------------------------------
#   * FreeSurfer 6.0.0 -> 8.2.0. VBG 2.0 requires FS >= 7.3: brain extraction
#     (-B 1) is mri_synthstrip, -P 4 is mri_synthseg, and the -M multi-atlas
#     block needs segment_subregions and mri_segment_hypothalamic_subunits.
#     None of those exist in FS 6.
#   * HD-BET removed entirely — VBG 2.0 Task 2.2 dropped both HD-BET code paths
#     in favour of SynthStrip (-B 1) / ANTs-BET (-B 2).
#   * Docker-in-Docker and the NVIDIA Container Toolkit removed. Both were
#     always questionable inside an image; under Apptainer they are unusable.
#   * Local CUDA toolkit (~4 GB) removed. PyTorch's own wheels ship the CUDA
#     runtime they need, and the host driver is injected at run time by
#     `docker --gpus all` / `apptainer --nv`. Nothing here compiles CUDA code.
#   * ANTs is built from a pinned commit instead of COPYed from the build host,
#     so the image is reproducible from the Dockerfile alone.
#   * MRtrix3 built from the CMake-era pin (the old `dev_sans_cmake` branch no
#     longer exists; upstream deleted ./configure && ./build in Oct 2023).
#   * FastSurfer installed via uv + .venv, NOT `conda env create -f
#     env/fastsurfer.yml`. That file is an empty placeholder in current
#     FastSurfer — the old recipe silently produced a broken install.
###############################################################################

ARG UBUNTU_TAG=24.04
# Pinned by DIGEST, not just tag. `ubuntu:24.04` is a moving tag that Canonical
# republishes for security updates, and BuildKit re-resolves it against the
# registry on every build rather than trusting whatever is in the local store.
# Confirmed live: a rebuild on 2026-09-05 resolved
# sha256:33ceb719... where the image it was replacing had been built from
# sha256:561618e2..., and that one difference invalidated EVERY downstream
# layer -- FSL, ANTs, MRtrix3 and FreeSurfer all rebuilt from source despite a
# 105 GB warm cache. The wasted hours are the small half of the problem: an
# unpinned base also means this image is NOT reproducible from this file alone,
# which the header above explicitly claims it is.
#
# The tag is kept next to the digest purely as documentation for a human
# reading `FROM`; the digest is what actually resolves.
#
# To adopt a newer base -- do it deliberately, e.g. for a CVE fix, and expect a
# full rebuild by design:
#   docker buildx imagetools inspect ubuntu:24.04 --format '{{.Manifest.Digest}}'
ARG UBUNTU_DIGEST=sha256:33ceb71981b602c1a7443a53469e4dba065f7503eab3078a2d7a57a2ab987517

###############################################################################
# Pinned versions — single source of truth, consumed by the stages below.
###############################################################################
# ANTs v2.4.4.post20.
ARG ANTS_COMMIT=40ee2d22

# MRtrix3: KUL fork of upstream dev. Carries one fix on top of dev (an
# f-string prefix missing from dwifslpreproc's mrconvert call). VBG itself
# only uses mrcalc/mrstats/labelconvert/mrhistmatch/mrthreshold/mrmath/
# mrfilter/maskfilter, none of which are affected — the fork is used here
# purely to keep the container identical to the validated bare-metal install.
ARG MRTRIX3_REPO=https://github.com/Rad-dude/mrtrix3.git
ARG MRTRIX3_BRANCH=fix/mrconvert-direct-io-arg-binding
ARG MRTRIX3_COMMIT=5a643594b0c0dae60f09e62cd6b1c56314a19354

# FreeSurfer 8.2.0. The ubuntu24 .deb matches the ubuntu:24.04 base below;
# if you change UBUNTU_TAG to 22.04, change this to freesurfer_ubuntu22-... too.
ARG FREESURFER_VERSION=8.2.0
ARG FREESURFER_DEB=freesurfer_ubuntu24-8.2.0_amd64.deb

# FSL 6.0.7.23 component versions come from the FSL conda channel; see the
# fsl stage for why this is a package subset rather than fslinstaller.py.
ARG FSL_FLAVOUR=subset

# FastSurfer. Pinned to the commit validated alongside FreeSurfer 8.2.0 here.
ARG FASTSURFER_COMMIT=700cfdc
ARG FASTSURFER_REPO=https://github.com/Deep-MI/FastSurfer.git

# PyTorch CUDA backend baked into FastSurfer's venv. This CANNOT be probed at
# build time (no GPU inside `docker build`), so it is a deliberate choice:
#   cu126  — default. Works on driver >= 525 via CUDA minor-version
#            compatibility, which covers essentially every current cluster.
#   cu128  — newer GPUs (Blackwell), needs a newer driver.
#   cu118  — very old drivers.
#   cpu    — no GPU anywhere, or you want a much smaller image.
# It must be a backend that publishes a wheel for the torch version pinned in
# FastSurfer's requirements.txt (2.7.1 at $FASTSURFER_COMMIT -> cu118/cu126/cu128).
ARG TORCH_BACKEND=cu126

# KUL_VBG itself. Commit-pinned for reproducibility; the branch is recorded so
# you can bump the commit to its tip without hunting for the branch name.
ARG VBG_REPO=https://github.com/KUL-Radneuron/KUL_VBG.git
ARG VBG_BRANCH=KUL_VBG_2.0
ARG VBG_COMMIT=56dc02a


###############################################################################
# Stage 1 — ANTs
#
# Built with gcc-12, NOT the ubuntu:24.04 default gcc-13. ANTs at this commit
# vendors a mid-2023 ITKv5 snapshot whose auto-generated Enum headers (e.g.
# itkMathematicalMorphologyEnums.h) rely on <cstdint> arriving transitively via
# another standard header. GCC 13 stopped doing that, so the ITK sub-build dies
# partway with "'uint8_t' was not declared". This is a documented ITK/GCC13
# incompatibility (ITK issue #4607), not something specific to this image.
# Building with gcc-12 sidesteps the entire class of breakage rather than
# patching vendored ITK sources one file at a time as each surfaces.
#
# Portability note: ANTs' own CMake sets -march=corei7 -mtune=native. Only
# -march restricts which CPUs the binary will RUN on, and corei7 is a
# Nehalem-era baseline that every current x86-64 server satisfies; -mtune is a
# scheduling hint that resolves to the build host's microarchitecture but does
# not emit instructions the target cannot execute. So an image built here is
# safe to ship to a different server -- which is the whole point of the SIF.
# Do not "optimise" this by adding -march=native: that WOULD bake in the build
# host's ISA and produce SIGILL on an older target.
###############################################################################
FROM ubuntu:${UBUNTU_TAG}@${UBUNTU_DIGEST} AS ants-builder
ARG ANTS_COMMIT
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        g++-12 \
        gcc-12 \
        git \
        make \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/ANTsX/ANTs.git /src/ANTs \
    && git -C /src/ANTs checkout "${ANTS_COMMIT}"

# The superbuild builds ITK first, then ANTS-build. `make install` must be run
# from the ANTS-build subdirectory — installing from the top-level build dir
# installs ITK instead of ANTs, which is the classic way to end up with an
# /opt/ANTs/bin that has no antsRegistration in it.
RUN mkdir -p /src/ANTs/build && cd /src/ANTs/build \
    && cmake \
        -DCMAKE_C_COMPILER=/usr/bin/gcc-12 \
        -DCMAKE_CXX_COMPILER=/usr/bin/g++-12 \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_TESTING=OFF \
        -DBUILD_SHARED_LIBS=OFF \
        -DCMAKE_INSTALL_PREFIX=/opt/ANTs \
        .. \
    && make -j"$(nproc)" \
    && cd ANTS-build && make install

# Fail the build here rather than at run time if the install landed wrong.
# These are exactly the ANTs entry points KUL_VBG.sh calls.
RUN set -eu; for _b in antsRegistration antsApplyTransforms Atropos ImageMath \
        N4BiasFieldCorrection ThresholdImage antsRegistrationSyN.sh \
        antsRegistrationSyNQuick.sh antsBrainExtraction.sh antsAtroposN4.sh; do \
        test -e "/opt/ANTs/bin/${_b}" || { echo "MISSING ANTs binary: ${_b}" >&2; exit 1; }; \
    done; echo "ANTs install verified"


###############################################################################
# Stage 2 — MRtrix3 (command-line only)
#
# -DMRTRIX_BUILD_GUI=OFF drops mrview and with it the whole Qt dependency tree;
# -DMRTRIX_ENABLE_GPU=OFF drops the Dawn/Slang shader stack that only the GUI
# uses. VBG calls no GUI tool (verified: only mrcalc, mrstats, labelconvert,
# mrhistmatch, mrthreshold, mrmath, mrfilter, maskfilter), and an image has no
# display anyway. This removes several GB and a long chunk of build time.
###############################################################################
FROM ubuntu:${UBUNTU_TAG}@${UBUNTU_DIGEST} AS mrtrix-builder
ARG MRTRIX3_REPO
ARG MRTRIX3_BRANCH
ARG MRTRIX3_COMMIT
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        git \
        libfftw3-dev \
        libpng-dev \
        libtiff5-dev \
        ninja-build \
        pkg-config \
        python3 \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone -b "${MRTRIX3_BRANCH}" "${MRTRIX3_REPO}" /src/mrtrix3 \
    && git -C /src/mrtrix3 checkout "${MRTRIX3_COMMIT}" \
    && git -C /src/mrtrix3 submodule update --init --recursive

# --prefix points back into the clone so bin/ and lib/mrtrix3/ land in the
# layout the rest of the KUL tooling expects (same as setup_environment.sh).
RUN cd /src/mrtrix3 \
    && cmake -B build -GNinja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/opt/mrtrix3 \
        -DMRTRIX_BUILD_GUI=OFF \
        -DMRTRIX_ENABLE_GPU=OFF \
        -DMRTRIX_BUILD_TESTS=OFF \
    && cmake --build build -j"$(nproc)" \
    && cmake --install build

RUN set -eu; for _b in mrcalc mrstats labelconvert mrhistmatch mrthreshold \
        mrmath mrfilter maskfilter mrconvert; do \
        test -x "/opt/mrtrix3/bin/${_b}" || { echo "MISSING MRtrix3 binary: ${_b}" >&2; exit 1; }; \
    done; echo "MRtrix3 install verified"


###############################################################################
# Stage 3 — FSL
#
# FSL 6.0.7 is a conda distribution, so this installs a dependency-resolved
# package SUBSET from FSL's own conda channel rather than running
# fslinstaller.py (which installs ~11 GB, most of it FSLeyes, FIRST models and
# atlases that VBG never touches).
#
# The subset is measured, not guessed. Every FSL call in the KUL_VBG tree:
#     fslmaths(114) fslstats(15) fslreorient2std(6) fslswapdim(2)
#     fslorient(2)                    -> fsl-avwutils
#     imcp(1)                         -> fslpy
#     MNI152_T1_1mm_brain.nii.gz (QC) -> fsl-data_standard
# VBG does brain extraction with SynthStrip/ANTs and tissue segmentation with
# ANTs Atropos, so it never calls FSL's own bet/fast/first.
#
# Because this goes through conda, the shared-library closure is resolved
# automatically — this is not a hand-copied pile of binaries, which is where
# subset installs normally break. If a rarely-hit path ever calls an FSL tool
# outside this set it fails loudly with "command not found" at that call,
# never silently or with wrong numbers.
#
# Build with --build-arg FSL_FLAVOUR=full to install the complete distribution
# instead; nothing else in this file needs to change.
###############################################################################
FROM ubuntu:${UBUNTU_TAG}@${UBUNTU_DIGEST} AS fsl-builder
ARG FSL_FLAVOUR
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        bzip2 ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*

# Micromamba is used purely as the solver/installer; it is not kept in the
# final image, and the resulting /opt/fsl tree is a plain conda prefix.
RUN curl -fsSL --retry 5 --retry-delay 5 --retry-all-errors \
        https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj -C /usr/local bin/micromamba

# The full-flavour package list is a manifest rather than an inline list: the
# FSL channel publishes no "fsl" metapackage, so "full" has to be enumerated,
# and enumerating it by hand invites silent drift. fsl-packages-full.txt is
# generated from a reference bare-metal install — see its header.
COPY fsl-packages-full.txt /tmp/fsl-packages-full.txt

RUN set -eux; \
    if [ "${FSL_FLAVOUR}" = "full" ]; then \
        _pkgs="$(sed -e 's/#.*//' -e '/^[[:space:]]*$/d' /tmp/fsl-packages-full.txt | tr '\n' ' ')"; \
    else \
        _pkgs="fsl-avwutils fslpy fsl-data_standard"; \
    fi; \
    echo "Installing FSL flavour '${FSL_FLAVOUR}': ${_pkgs}"; \
    micromamba create -y -p /opt/fsl \
        -c https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/public \
        -c conda-forge \
        python=3.12 ${_pkgs}; \
    micromamba clean -y --all; \
    rm -f /tmp/fsl-packages-full.txt

RUN set -eu; for _b in fslmaths fslstats fslreorient2std fslswapdim fslorient imcp; do \
        test -x "/opt/fsl/bin/${_b}" || { echo "MISSING FSL binary: ${_b}" >&2; exit 1; }; \
    done; \
    test -f /opt/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz \
        || { echo "MISSING FSL standard-space MNI152 template" >&2; exit 1; }; \
    echo "FSL install verified"


###############################################################################
# Stage 4 — runtime
###############################################################################
FROM ubuntu:${UBUNTU_TAG}@${UBUNTU_DIGEST}

LABEL maintainer="radwanphd@gmail.com" \
      org.opencontainers.image.title="KUL_VBG" \
      org.opencontainers.image.description="Virtual Brain Grafting — whole-brain parcellation in the presence of large lesions" \
      org.opencontainers.image.source="https://github.com/KUL-Radneuron/KUL_VBG" \
      org.opencontainers.image.licenses="MPL-2.0"

ARG FREESURFER_VERSION
ARG FREESURFER_DEB
ARG FASTSURFER_REPO
ARG FASTSURFER_COMMIT
ARG TORCH_BACKEND
ARG VBG_REPO
ARG VBG_BRANCH
ARG VBG_COMMIT

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

# Runtime dependencies only — no compilers, no Qt, no CUDA toolkit.
# tcsh/bc/perl are FreeSurfer's; libgomp1 is used by both ANTs and FSL.
RUN apt-get update && apt-get install -y --no-install-recommends \
        bc \
        ca-certificates \
        curl \
        file \
        git \
        libfftw3-single3 \
        libfftw3-double3 \
        libgomp1 \
        libpng16-16t64 \
        libquadmath0 \
        libtiff6 \
        locales \
        perl \
        python3 \
        python3-venv \
        tcsh \
        tzdata \
        unzip \
        wget \
        zlib1g \
    && locale-gen en_US.UTF-8 \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# FreeSurfer 8.2.0
#
# Installed as a real .deb via dpkg + `apt-get install -f`, NOT unpacked with
# `dpkg-deb -x`. A payload-only extraction skips the .deb's declared runtime
# library dependencies, which breaks freeview and other ITK-linked tools
# ("libITKCommon-5.3.so.1: cannot open shared object file") — ITK is not
# bundled inside FreeSurfer's tree, only referenced as a dependency apt is
# expected to supply. `dpkg -i` exiting non-zero here is expected (it reports
# the missing deps but still unpacks); `apt-get install -f` then resolves them
# and finishes configuring the package.
#
# NO LICENSE FILE IS BAKED IN. FreeSurfer's licence terms are per-user, so the
# licence must be supplied at run time as a bind mount — see entrypoint.sh.
# ---------------------------------------------------------------------------
RUN set -eux; \
    curl -fSL --retry 5 --retry-delay 5 --retry-all-errors \
        -o /tmp/freesurfer.deb \
        "https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/${FREESURFER_VERSION}/${FREESURFER_DEB}"; \
    apt-get update; \
    dpkg -i /tmp/freesurfer.deb || true; \
    apt-get install -f -y --no-install-recommends; \
    rm -f /tmp/freesurfer.deb; \
    rm -rf /var/lib/apt/lists/*

ENV FREESURFER_HOME=/usr/local/freesurfer/${FREESURFER_VERSION}
# Same value under a name the host environment will never define. Apptainer
# merges the host environment into the container by default, so a user whose
# login shell exports its own FREESURFER_HOME could otherwise shadow the
# container's. entrypoint.sh reads this one, which keeps it both host-proof and
# automatically correct if FREESURFER_VERSION is ever bumped.
ENV KUL_VBG_FREESURFER_HOME=${FREESURFER_HOME}
ENV PATH="${FREESURFER_HOME}/bin:${PATH}" \
    SUBJECTS_DIR=/data \
    FS_LICENSE=${FREESURFER_HOME}/license.txt \
    FSFAST_HOME=${FREESURFER_HOME}/fsfast \
    MNI_DIR=${FREESURFER_HOME}/mni

# ---------------------------------------------------------------------------
# FreeSurfer 8.2.0 official post-release patches — REQUIRED, not optional.
#
# The 8.2.0 .deb ships samseg's compiled gems bindings built ONLY for Python
# 3.12 (gemsbindings.cpython-312-*.so), while FreeSurfer's own bundled
# interpreter is Python 3.8.13. Python 3.8 cannot load a cpython-312 extension,
# so `segment_subregions` dies on import with
#     ModuleNotFoundError: No module named 'samseg.gems.gemsbindings'
# and every subregion segmentation silently fails: hippocampal subfields,
# amygdala nuclei, thalamic nuclei and brainstem substructures. Confirmed on a
# real run — VBG logs a non-fatal WARNING and skips them, so an unpatched image
# looks like it succeeded while quietly producing four fewer segmentations.
#
# Upstream ships the fix through fs820_updates.sh; the patch's own note for
# this file reads "add gemsbindings.cpython python 3.8 library". The ubuntu24
# patch set also replaces recon-all, mri_glmfit, rca-surfreg, segmentHA_T2.sh
# and others, so this is not a one-file workaround — an unpatched install is
# simply an incomplete FreeSurfer.
#
# The archive is PINNED rather than taking whatever "default" currently points
# at, so a rebuild reproduces this exact FreeSurfer. Bump FREESURFER_PATCH to
# adopt a newer patch set deliberately. (The script still fetches "default"
# first, so the build needs network either way; passing the archive as $1 makes
# it use the pinned one regardless of what default says.)
#
# wget, not curl: this host 403s curl's default user-agent but serves wget fine.
#
# `yes | tr -d '\n'` is not decoration — the script prompts per file with
#     read -p "Update ...? [Yy/Nn]: " -n 1 -r
# and there is no batch/force flag. Two traps:
#   * With no stdin (a plain docker build) every read gets EOF, REPLY stays
#     empty, and the script cheerfully reports "Skipping replacing ..." for all
#     11 files while still exiting 0. The patch appears to run and does nothing.
#   * Plain `yes |` is WORSE THAN IT LOOKS: yes emits "y\n", `read -n 1` takes
#     the "y", then the NEXT read consumes the leftover "\n", which fails the
#     ^[Yy]$ test. Files then alternate applied/skipped — observed live, 5 of 11
#     applied and the gems binding among those skipped.
# Stripping newlines yields an unbroken "yyyy..." stream so every prompt gets a
# real "y". Verified: 11 applied, 0 skipped, segment_subregions runs.
#
# The verification step below is what catches a regression here — it asserts the
# .so exists AND that segment_subregions actually executes.
# ---------------------------------------------------------------------------
ARG FREESURFER_PATCH=patch_data_20260625.tgz
RUN set -eux; \
    _pd=/tmp/fspatch; mkdir -p "${_pd}"; cd "${_pd}"; \
    wget -q "https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/${FREESURFER_VERSION}/fs820_updates.sh"; \
    chmod +x fs820_updates.sh; \
    yes | tr -d '\n' | ./fs820_updates.sh "${FREESURFER_PATCH}"; \
    cd /; rm -rf "${_pd}"

# Confirm the .deb landed where the ENV above says, and that the FS 7.3+ tools
# VBG depends on are actually present. Cheap here, painful to discover mid-run.
RUN set -eu; \
    test -f "${FREESURFER_HOME}/SetUpFreeSurfer.sh" \
        || { echo "FreeSurfer not at ${FREESURFER_HOME} — check FREESURFER_DEB/FREESURFER_VERSION" >&2; \
             ls -d /usr/local/freesurfer* 2>/dev/null >&2; exit 1; }; \
    for _b in recon-all mri_synthstrip mri_synthseg segment_subregions \
              mri_segment_hypothalamic_subunits mri_surf2surf mri_aparc2aseg \
              mris_register mri_convert mri_mask; do \
        test -x "${FREESURFER_HOME}/bin/${_b}" \
            || { echo "MISSING FreeSurfer binary: ${_b}" >&2; exit 1; }; \
    done; \
    test -d "${FREESURFER_HOME}/subjects/fsaverage" \
        || { echo "MISSING fsaverage (needed by the -M Lausanne/Glasser block)" >&2; exit 1; }; \
    _gems="${FREESURFER_HOME}/python/lib/python3.8/site-packages/samseg/gems"; \
    test -f "${_gems}/gemsbindings.cpython-38-x86_64-linux-gnu.so" \
        || { echo "MISSING samseg cpython-38 gems bindings — the fs820 patch did not apply." >&2; \
             echo "  Without it segment_subregions fails and the -M block silently produces no" >&2; \
             echo "  hippocampal/amygdala/thalamic/brainstem segmentations." >&2; exit 1; }; \
    "${FREESURFER_HOME}/bin/segment_subregions" --help >/dev/null 2>&1 \
        || { echo "segment_subregions is present but does not run (samseg import failure?)" >&2; exit 1; }; \
    echo "FreeSurfer ${FREESURFER_VERSION} verified (incl. fs820 patches + segment_subregions)"

# ---------------------------------------------------------------------------
# ANTs, MRtrix3, FSL — from the builder stages
# ---------------------------------------------------------------------------
COPY --from=ants-builder   /opt/ANTs    /opt/ANTs
COPY --from=mrtrix-builder /opt/mrtrix3 /opt/mrtrix3
COPY --from=fsl-builder    /opt/fsl     /opt/fsl

ENV ANTSPATH=/opt/ANTs/bin \
    FSLDIR=/opt/fsl \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    FSLTCLSH=/opt/fsl/bin/fsltclsh \
    FSLWISH=/opt/fsl/bin/fslwish
# $FSLDIR/bin holds the real binaries (the conda packages install there);
# $FSLDIR/share/fsl/bin holds absolute-path wrapper shims generated by
# fsl-base's post-link script. bin/ is listed first so the image does not
# depend on micromamba having run those post-link scripts.
ENV PATH="/opt/ANTs/bin:/opt/mrtrix3/bin:/opt/fsl/bin:/opt/fsl/share/fsl/bin:${PATH}"

# ---------------------------------------------------------------------------
# FastSurfer (for KUL_VBG.sh -P 2 and -P 3)
#
# Installed with uv into a per-checkout .venv, which is FastSurfer's own
# supported native-install route now — env/fastsurfer.yml is an empty
# placeholder in current FastSurfer and `conda env create` against it produces
# a broken install.
#
# The `uv pip compile | uv pip sync` pipeline is deliberate and must not be
# collapsed into a bare `uv pip sync requirements.txt`: requirements.txt
# excludes the nvidia-* CUDA runtime packages from its pinned list (per its own
# header comment) and expects the compile step to pull them in transitively.
# `sync` alone installs only what is literally listed, silently omitting the
# runtime libs and leaving torch unable to import at run time
# ("libcudart.so.*: cannot open shared object file").
# ---------------------------------------------------------------------------
ENV FASTSURFER_HOME=/opt/FastSurfer
ENV PATH="${FASTSURFER_HOME}:${PATH}"

RUN curl -fsSL --retry 5 --retry-delay 5 --retry-all-errors \
        https://astral.sh/uv/install.sh | env UV_INSTALL_DIR=/usr/local/bin sh

RUN set -eux; \
    git clone "${FASTSURFER_REPO}" "${FASTSURFER_HOME}"; \
    git -C "${FASTSURFER_HOME}" checkout "${FASTSURFER_COMMIT}"; \
    cd "${FASTSURFER_HOME}"; \
    uv venv --python 3.12; \
    uv pip compile --torch-backend "${TORCH_BACKEND}" requirements.txt \
        --python .venv/bin/python \
    | uv pip sync --torch-backend "${TORCH_BACKEND}" --python .venv/bin/python -; \
    # uv keeps every downloaded wheel under /root/.cache/uv — ~5.8 GB here,
    # dominated by torch and the nvidia-* CUDA runtime packages. The .venv
    # already holds the installed copies, so the cache is pure dead weight in
    # every layer below this one. Must be cleaned in THIS RUN: a later
    # `rm -rf` would only mask it, leaving the bytes in this layer and in the
    # SIF. FastSurfer's own .git is deliberately KEPT (94 MB) — version.py
    # probes for it to report the checkout version, and there is no VERSION
    # file fallback, so dropping it would leave provenance unknown.
    uv cache clean; \
    rm -rf /root/.cache/uv

# Wrapper shim: run_fastsurfer.sh calls plain python3 and never activates its
# own venv, and KUL_VBG.sh invokes it by bare name via PATH (`which
# run_fastsurfer.sh`). Renaming the real script and dropping a wrapper at the
# original filename covers both invocation styles without patching VBG source.
RUN set -eux; \
    mv "${FASTSURFER_HOME}/run_fastsurfer.sh" "${FASTSURFER_HOME}/run_fastsurfer.real.sh"; \
    printf '%s\n' \
        '#!/usr/bin/env bash' \
        '# Wrapper shim: activates FastSurfer'"'"'s uv-managed .venv, then hands off.' \
        '# Scoped to this process via exec — never leaks back to the calling shell.' \
        'set -euo pipefail' \
        'FASTSURFER_HOME="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" &>/dev/null && pwd)"' \
        'source "$FASTSURFER_HOME/.venv/bin/activate"' \
        'exec "$FASTSURFER_HOME/run_fastsurfer.real.sh" "$@"' \
        > "${FASTSURFER_HOME}/run_fastsurfer.sh"; \
    chmod +x "${FASTSURFER_HOME}/run_fastsurfer.sh"

# Bake the network checkpoints into the image. FastSurfer otherwise fetches
# them on first run, which fails on a compute node with no outbound network
# and/or tries to write into a read-only SIF. This is the single most important
# difference between an image that works on an HPC node and one that does not.
#
# PYTHONPATH is REQUIRED and easy to lose: download_checkpoints.py does an
# absolute `from FastSurferCNN.utils.checkpoint import ...`, but running it as
# `python FastSurferCNN/download_checkpoints.py` puts the SCRIPT's directory
# ($FASTSURFER_HOME/FastSurferCNN) on sys.path, not $FASTSURFER_HOME — so the
# package cannot be found and the step dies with ModuleNotFoundError. Nothing
# pip-installs FastSurfer into its own venv, so there is no other route.
# (run_fastsurfer.real.sh sets PYTHONPATH itself, so only this build-time
# invocation needs it — the runtime path was never affected.)
#
# Beware when testing this by hand: a PYTHONPATH with a trailing ':' contains an
# empty entry, which Python reads as "current directory" — so it appears to work
# from inside $FASTSURFER_HOME and hides the bug. Test with `env -u PYTHONPATH`.
RUN set -eux; \
    cd "${FASTSURFER_HOME}"; \
    PYTHONPATH="${FASTSURFER_HOME}" \
        .venv/bin/python FastSurferCNN/download_checkpoints.py --all; \
    _n="$(find "${FASTSURFER_HOME}/checkpoints" -name '*.pkl' 2>/dev/null | wc -l)"; \
    test "${_n}" -gt 0 \
        || { echo "download_checkpoints.py exited 0 but wrote no .pkl files" >&2; exit 1; }; \
    echo "FastSurfer checkpoints baked in: ${_n} files"

# ---------------------------------------------------------------------------
# KUL_VBG
# ---------------------------------------------------------------------------
# share/Test_data (~105 MB of defaced example T1s) is deliberately KEPT: it is
# already public in the repo, it is negligible against the image size, and it
# lets the container be verified end-to-end without the user supplying data.
# Only .git is dropped.
RUN set -eux; \
    git clone -b "${VBG_BRANCH}" "${VBG_REPO}" /opt/KUL_VBG; \
    git -C /opt/KUL_VBG checkout "${VBG_COMMIT}"; \
    rm -rf /opt/KUL_VBG/.git; \
    chmod +x /opt/KUL_VBG/*.sh /opt/KUL_VBG/*.py

ENV PATH="/opt/KUL_VBG:${PATH}"

# Python packages for KUL_VBG_QC.py, KUL_lesion_overlap.py and
# remap_lausanne_to_msbp.py. Deliberately a plain venv rather than another
# conda env: these are four pure-wheel packages and this keeps the layer small.
ENV VBG_PYENV=/opt/vbg-pyenv
# /usr/bin/python3 explicitly, NOT bare `python3`: by this point /opt/fsl/bin is
# already ahead of /usr/bin on PATH, so a bare `python3` resolves to FSL's conda
# interpreter and the venv silently gets built on top of it. That does work, but
# it couples the QC tooling to whatever Python the FSL packages happen to pin —
# change FSL_FLAVOUR or bump the FSL pin and this venv's base can move under it.
RUN /usr/bin/python3 -m venv "${VBG_PYENV}" \
    && "${VBG_PYENV}/bin/pip" install --no-cache-dir --upgrade pip \
    && "${VBG_PYENV}/bin/pip" install --no-cache-dir \
        "numpy" "scipy" "nibabel" "matplotlib"
ENV PATH="${VBG_PYENV}/bin:${PATH}"

# ---------------------------------------------------------------------------
# Runtime environment for non-interactive shells
#
# KUL_VBG.sh is normally invoked directly (not through a login shell), so the
# ENVs above do the real work. These profile.d files exist so that
# `docker run ... /bin/bash` and `apptainer shell` also land in a working
# environment, and so FreeSurfer's own SetUpFreeSurfer.sh gets sourced for the
# handful of variables it sets beyond FREESURFER_HOME.
# ---------------------------------------------------------------------------
RUN set -eux; \
    printf '%s\n' \
        'export FREESURFER_HOME='"${FREESURFER_HOME}" \
        'export SUBJECTS_DIR=${SUBJECTS_DIR:-/data}' \
        'source $FREESURFER_HOME/SetUpFreeSurfer.sh >/dev/null 2>&1 || true' \
        > /etc/profile.d/z10-freesurfer.sh; \
    printf '%s\n' \
        'export FSLDIR=/opt/fsl' \
        'export FSLOUTPUTTYPE=NIFTI_GZ' \
        'export PATH=$FSLDIR/bin:$FSLDIR/share/fsl/bin:$PATH' \
        > /etc/profile.d/z20-fsl.sh; \
    printf '%s\n' \
        'export ANTSPATH=/opt/ANTs/bin' \
        'export PATH=$ANTSPATH:$PATH' \
        > /etc/profile.d/z30-ants.sh; \
    printf '%s\n' \
        'export PATH=/opt/mrtrix3/bin:$PATH' \
        > /etc/profile.d/z40-mrtrix3.sh; \
    printf '%s\n' \
        'export FASTSURFER_HOME=/opt/FastSurfer' \
        'export PATH=/opt/KUL_VBG:/opt/vbg-pyenv/bin:$FASTSURFER_HOME:$PATH' \
        > /etc/profile.d/z50-kul_vbg.sh

# ---------------------------------------------------------------------------
# FreeSurfer licence plumbing
#
# The licence is NEVER baked into the image (FreeSurfer's terms are per-user);
# it is supplied at run time. Getting this right under Apptainer needs care:
# a SIF is READ-ONLY, so the old entrypoint's `cp <mounted licence>
# $FREESURFER_HOME/license.txt` cannot work there.
#
# Instead, $FREESURFER_HOME/license.txt is baked in as a SYMLINK to the stable
# path /licence/license.txt, and users bind their licence onto that. Nothing is
# ever written inside the image:
#   - KUL_VBG.sh's own check (FS_lic_search[0] = $FREESURFER_HOME/license.txt)
#     follows the symlink and passes, since [[ -f ]] resolves symlinks.
#   - FreeSurfer's and FastSurfer's internal lookups resolve it the same way.
#   - entrypoint.sh additionally exports FS_LICENSE, which is the last entry in
#     that same search array, so an alternative bind path also works.
#
# The placeholder file exists so that the bind target is present: binding a
# file onto a path that does not exist inside a read-only SIF requires
# overlay/underlay support that is not guaranteed on every cluster.
# entrypoint.sh rejects the placeholder by checking the licence is non-empty,
# so an unbound licence fails with a clear message instead of a cryptic
# FreeSurfer error twenty minutes into a run.
# ---------------------------------------------------------------------------
RUN mkdir -p /licence \
    && : > /licence/license.txt \
    && chmod 0644 /licence/license.txt \
    && ln -sf /licence/license.txt "${FREESURFER_HOME}/license.txt"

# /data is the conventional bind point for input/output. Created here so it
# exists inside a read-only SIF even when nothing is bound over it.
RUN mkdir -p /data
WORKDIR /data

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
CMD ["KUL_VBG.sh", "-h"]
