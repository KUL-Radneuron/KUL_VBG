#!/usr/bin/env bash
###############################################################################
# KUL_VBG 2.0 — container build helper
#
# Builds the Docker image and, optionally, converts it to an Apptainer SIF for
# use on a rootless multi-user server.
#
# Usage:
#   ./build.sh                       # Docker image only
#   ./build.sh --sif                 # Docker image, then convert to SIF
#   ./build.sh --sif-only            # convert an already-built image to SIF
#   ./build.sh --torch-backend cu128 # GPU backend for FastSurfer's torch
#   ./build.sh --fsl full            # full FSL instead of the VBG subset
#   ./build.sh --push                # push to the registry after building
#
# Run ./build.sh --help for the full list.
###############################################################################
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

# ── Defaults ─────────────────────────────────────────────────────────────────
IMAGE_NAME="kul_vbg"
IMAGE_TAG="2.0"
REGISTRY="radwankul"          # for --push
TORCH_BACKEND="cu126"
FSL_FLAVOUR="subset"
SIF_PATH=""                   # defaults to ./KUL_VBG_<tag>.sif
DO_DOCKER=1
DO_SIF=0
DO_PUSH=0
NO_CACHE=0
EXTRA_BUILD_ARGS=()

# ── Colours (only when attached to a terminal) ───────────────────────────────
if [ -t 1 ]; then
    C_OK=$'\033[32m'; C_WARN=$'\033[33m'; C_ERR=$'\033[31m'
    C_HEAD=$'\033[1;36m'; C_OFF=$'\033[0m'
else
    C_OK=""; C_WARN=""; C_ERR=""; C_HEAD=""; C_OFF=""
fi
log()  { echo "${C_HEAD}==>${C_OFF} $*"; }
ok()   { echo "${C_OK}[ok]${C_OFF} $*"; }
warn() { echo "${C_WARN}[warn]${C_OFF} $*" >&2; }
die()  { echo "${C_ERR}[error]${C_OFF} $*" >&2; exit 1; }

usage() {
    sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    cat <<EOF

Options:
  --name NAME            image name             (default: ${IMAGE_NAME})
  --tag TAG              image tag              (default: ${IMAGE_TAG})
  --registry NS          registry namespace for --push (default: ${REGISTRY})
  --torch-backend TAG    cu118|cu126|cu128|cpu  (default: ${TORCH_BACKEND})
  --fsl subset|full      FSL install flavour    (default: ${FSL_FLAVOUR})
  --sif                  also build the Apptainer SIF
  --sif-only             skip the Docker build, only convert to SIF
  --sif-path PATH        output SIF path        (default: ./KUL_VBG_<tag>.sif)
  --push                 docker push after a successful build
  --no-cache             docker build --no-cache
  --build-arg K=V        extra --build-arg, repeatable
  -h, --help             this message

Notes on --torch-backend: this cannot be detected automatically, because
'docker build' has no GPU. cu126 (the default) covers driver >= 525 through
CUDA minor-version compatibility, which is almost every current cluster. Use
cu128 for Blackwell-era GPUs on a new driver, cu118 for very old drivers, or
cpu for a much smaller image with no GPU support. Check the target server with
'nvidia-smi' -- the "CUDA Version" it reports is the highest your driver
supports.
EOF
}

# ── Arguments ────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --name)          IMAGE_NAME="$2"; shift 2 ;;
        --tag)           IMAGE_TAG="$2"; shift 2 ;;
        --registry)      REGISTRY="$2"; shift 2 ;;
        --torch-backend) TORCH_BACKEND="$2"; shift 2 ;;
        --fsl)           FSL_FLAVOUR="$2"; shift 2 ;;
        --sif)           DO_SIF=1; shift ;;
        --sif-only)      DO_SIF=1; DO_DOCKER=0; shift ;;
        --sif-path)      SIF_PATH="$2"; shift 2 ;;
        --push)          DO_PUSH=1; shift ;;
        --no-cache)      NO_CACHE=1; shift ;;
        --build-arg)     EXTRA_BUILD_ARGS+=(--build-arg "$2"); shift 2 ;;
        -h|--help)       usage; exit 0 ;;
        *)               die "unknown option: $1  (try --help)" ;;
    esac
done

IMAGE_REF="${IMAGE_NAME}:${IMAGE_TAG}"
: "${SIF_PATH:=KUL_VBG_${IMAGE_TAG}.sif}"

case "${FSL_FLAVOUR}" in
    subset|full) ;;
    *) die "--fsl must be 'subset' or 'full', got '${FSL_FLAVOUR}'" ;;
esac

case "${TORCH_BACKEND}" in
    cu118|cu126|cu128|cpu) ;;
    *) warn "unusual --torch-backend '${TORCH_BACKEND}'. It must be a backend that
       publishes a wheel for the torch version FastSurfer pins, or the build
       will fail during 'uv pip sync'." ;;
esac

# ── Docker build ─────────────────────────────────────────────────────────────
if [[ ${DO_DOCKER} -eq 1 ]]; then
    command -v docker >/dev/null 2>&1 || die "docker not found on PATH"
    docker info >/dev/null 2>&1 || die "cannot talk to the Docker daemon (is it running, and are you in the 'docker' group?)"

    log "Building ${IMAGE_REF}"
    echo "      FSL flavour   : ${FSL_FLAVOUR}"
    echo "      torch backend : ${TORCH_BACKEND}"

    # BuildKit matters here, and not just for speed. The ants/mrtrix/fsl stages
    # are independent, so BuildKit builds them in PARALLEL and skips any stage
    # the target does not need. The legacy builder does neither: it walks every
    # stage sequentially, so even `--target fsl-builder` compiles ANTs first.
    # The build still succeeds without it, just markedly slower.
    if docker buildx version >/dev/null 2>&1; then
        BUILDER="buildx"
        echo "      builder       : BuildKit (buildx)"
    else
        BUILDER="legacy"
        echo "      builder       : legacy (deprecated)"
    fi
    echo

    if [[ "${BUILDER}" = "legacy" ]]; then
        warn "docker buildx is not installed, so the deprecated legacy builder will be"
        warn "used. The build works, but stages run sequentially instead of in parallel."
        warn "To speed it up considerably:  sudo apt-get install -y docker-buildx"
        echo
    fi

    warn "This builds ANTs and MRtrix3 from source and downloads FreeSurfer 8.2.0"
    warn "plus FastSurfer's checkpoints. Expect 1-3 hours and ~40 GB of scratch"
    warn "space on the first build; later builds reuse the layer cache."
    echo

    build_args=(
        --build-arg "FSL_FLAVOUR=${FSL_FLAVOUR}"
        --build-arg "TORCH_BACKEND=${TORCH_BACKEND}"
    )

    # REQUIRED for the Apptainer conversion, not an optimisation.
    #
    # By default BuildKit attaches provenance/SBOM attestations, which forces
    # the result into a manifest LIST containing (a) the real image manifest and
    # (b) an attestation manifest whose "layer" is an in-toto JSON blob, not a
    # tar. `apptainer build docker-daemon://` walks that list and tries to
    # untar the JSON, failing with a thoroughly misleading
    #     archive/tar: invalid tar header
    # that looks like a corrupt image but is nothing of the sort. Disabling both
    # attestations makes BuildKit export a single plain manifest again.
    #
    # These flags exist only under buildx; the legacy builder never produced
    # attestations in the first place, so it needs (and rejects) them.
    if [[ "${BUILDER}" = "buildx" ]]; then
        build_args+=(--provenance=false --sbom=false)
    fi
    [[ ${NO_CACHE} -eq 1 ]] && build_args+=(--no-cache)
    [[ ${#EXTRA_BUILD_ARGS[@]} -gt 0 ]] && build_args+=("${EXTRA_BUILD_ARGS[@]}")

    docker build "${build_args[@]}" -t "${IMAGE_REF}" -f Dockerfile . \
        || die "docker build failed"

    ok "built ${IMAGE_REF}  ($(docker image inspect "${IMAGE_REF}" --format '{{.Size}}' \
        | awk '{printf "%.1f GB", $1/1024/1024/1024}'))"

    # Smoke test. `KUL_VBG.sh -h` exits non-zero by design (it is the usage
    # path), and the entrypoint refuses to run without a licence, so neither is
    # a useful check. Verifying the tool inventory instead confirms every
    # dependency actually landed on PATH inside the final image.
    log "Smoke-testing the image"
    docker run --rm --entrypoint /bin/bash "${IMAGE_REF}" -lc '
        set -e
        for t in recon-all mri_synthstrip mri_synthseg segment_subregions \
                 antsRegistration Atropos ImageMath antsRegistrationSyN.sh \
                 mrcalc mrstats labelconvert mrhistmatch \
                 fslmaths fslstats fslreorient2std imcp \
                 run_fastsurfer.sh KUL_VBG.sh KUL_VBG_multiparc.sh; do
            command -v "$t" >/dev/null || { echo "MISSING: $t"; exit 1; }
        done
        python3 -c "import numpy, scipy, nibabel, matplotlib" \
            || { echo "MISSING: python QC dependencies"; exit 1; }
        test -s /opt/KUL_VBG/atlasses/New/atlases/glasser/lh.HCPMMP1.annot \
            || { echo "MISSING: Glasser annot data"; exit 1; }
        # Must actually RUN, not merely exist. The unpatched FreeSurfer 8.2.0
        # .deb ships this binary but only the cpython-312 samseg bindings, so it
        # dies on import under the bundled Python 3.8. A command -v check passes
        # while every subregion segmentation (hippocampal, amygdala, thalamic,
        # brainstem) silently fails at run time. That is how a broken image
        # shipped once. NOTE: no apostrophes in here, this whole block is inside
        # a single-quoted bash -lc string.
        segment_subregions --help >/dev/null 2>&1 \
            || { echo "BROKEN: segment_subregions does not run (fs820 patch missing?)"; exit 1; }
        # Checkpoints must be baked in: fetching them at run time fails on an
        # offline compute node and cannot write into a read-only SIF.
        test -n "$(find /opt/FastSurfer/checkpoints -name "*.pkl" 2>/dev/null | head -1)" \
            || { echo "MISSING: FastSurfer checkpoints"; exit 1; }
        # torch must actually import — a wrong --torch-backend leaves the
        # nvidia-* runtime libs missing and only shows up here.
        /opt/FastSurfer/.venv/bin/python -c "import torch" \
            || { echo "MISSING: FastSurfer torch runtime"; exit 1; }
        echo "all tools present"
    ' || die "smoke test failed — the image is not usable"
    ok "smoke test passed"
fi

# ── Apptainer conversion ─────────────────────────────────────────────────────
if [[ ${DO_SIF} -eq 1 ]]; then
    APPTAINER_BIN="$(command -v apptainer || command -v singularity || true)"
    [[ -n "${APPTAINER_BIN}" ]] || die "neither apptainer nor singularity found on PATH"

    log "Converting ${IMAGE_REF} -> ${SIF_PATH}"
    warn "This unpacks and re-compresses the whole image; it takes a while and"
    warn "needs roughly 2x the image size in \$APPTAINER_TMPDIR (currently: ${APPTAINER_TMPDIR:-/tmp})."

    # Why this goes via a verified `docker save` archive instead of straight
    # from docker-daemon://
    #
    # Docker's layer export was observed (2026-08-14, docker 29.1.3) to
    # INTERMITTENTLY emit a corrupt blob for this image's ~11 GB FreeSurfer
    # layer: correct length, wrong bytes, so the blob's content no longer
    # hashed to the digest it was stored under. Confirmed independently with
    # GNU tar ("gzip: invalid compressed data--crc error" partway through) and
    # by hashing the blob. Re-exporting the very same image produced a blob
    # that verified clean, so the stored layer was never damaged — only the
    # export was, and only sometimes. Smaller blobs in the same export were
    # fine every time.
    #
    # apptainer surfaces this as "archive/tar: invalid tar header", which reads
    # like a broken image and is thoroughly misleading — and you only find out
    # ~20 minutes into the conversion. Exporting to an archive first lets the
    # blobs be checked against their own digests in a couple of minutes, and a
    # bad export simply gets retried. Content-addressed blobs make this a free,
    # exact check: no reference copy needed, the name IS the checksum.
    _archive="${SIF_PATH%.sif}.docker.tar"
    _attempt=1
    _max_attempts=3
    while : ; do
        log "Exporting image to ${_archive} (attempt ${_attempt}/${_max_attempts})"
        docker save "${IMAGE_REF}" -o "${_archive}" || die "docker save failed"

        log "Verifying exported layer digests"
        _bad=0
        while read -r _blob; do
            [[ -n "${_blob}" ]] || continue
            _actual="$(tar -xOf "${_archive}" "blobs/sha256/${_blob}" 2>/dev/null \
                       | sha256sum | cut -d' ' -f1)"
            if [[ "${_actual}" != "${_blob}" ]]; then
                warn "corrupt layer: ${_blob:0:16}... hashed to ${_actual:0:16}..."
                _bad=$((_bad + 1))
            fi
        done < <(tar -tf "${_archive}" 2>/dev/null | grep '^blobs/sha256/' | sed 's|blobs/sha256/||')

        if [[ ${_bad} -eq 0 ]]; then
            ok "all exported layers verified"
            break
        fi

        rm -f "${_archive}"
        if [[ ${_attempt} -ge ${_max_attempts} ]]; then
            die "docker save produced corrupt layers ${_max_attempts} times running.

That is more than bad luck. Two things worth checking:
  * RAM. Intermittent corruption of only the largest object, with no kernel I/O
    errors logged, is a classic non-ECC bit-flip signature. Run memtest86+.
  * Disk space and the filesystem under \$TMPDIR and /var/lib/docker."
        fi
        warn "retrying the export"
        _attempt=$((_attempt + 1))
    done

    if "${APPTAINER_BIN}" build --force "${SIF_PATH}" "docker-archive://${_archive}"; then
        rm -f "${_archive}"
    else
        rm -f "${_archive}"
        warn "apptainer could not unpack the OCI layers — falling back to the flatten route."
        echo

        # Why a fallback exists at all.
        #
        # apptainer 1.5.3 cannot unpack this image's ~11 GB FreeSurfer layer. It
        # dies at tar entry 33,830 of 107,050 with "archive/tar: invalid tar
        # header", on an ordinary 1.1 MB file. The layer is provably fine: its
        # blob hashes to exactly its own digest, and GNU tar reads all 107,050
        # entries with exit 0. So this is an apptainer-side defect, not a bad
        # image — which is why the verified-archive path above cannot fix it.
        #
        # The way around it is to never hand apptainer a layer: `docker export`
        # flattens the image to a single uncompressed tar, GNU tar extracts it,
        # and apptainer only squashes a directory.
        #
        # The fast path is still tried first, so this self-heals if a future
        # apptainer fixes the bug.
        _work="$(dirname "${SIF_PATH}")/.vbg_sifbuild"
        _rootfs="${_work}/rootfs.tar"
        _sandbox="${_work}/sandbox"
        _cname="vbg_export_$$"

        warn "This needs roughly 2x the image size again in $(dirname "${SIF_PATH}")."
        rm -rf "${_work}"; mkdir -p "${_sandbox}"

        log "Flattening ${IMAGE_REF} with docker export"
        docker rm -f "${_cname}" >/dev/null 2>&1 || true
        docker create --name "${_cname}" "${IMAGE_REF}" >/dev/null || die "docker create failed"
        docker export "${_cname}" -o "${_rootfs}" || { docker rm -f "${_cname}" >/dev/null 2>&1; die "docker export failed"; }
        docker rm -f "${_cname}" >/dev/null 2>&1 || true

        log "Extracting the rootfs with GNU tar"
        # `unshare -r` alone maps only ONE uid to root, so tar aborts the moment
        # it meets any other uid in the image (6, 12, 43 ...). --map-auto maps
        # the whole /etc/subuid range. Without a subuid range, extract as the
        # current user: mode bits are preserved either way and every binary here
        # is world-readable, so only the recorded owner differs.
        # tar's exit status is deliberately tolerated — it returns non-zero for
        # ownership warnings on paths like var/log and var/cache/man that
        # nothing in this pipeline uses. The rootfs check below is the real gate.
        if unshare --map-auto --map-root-user true 2>/dev/null; then
            unshare --map-auto --map-root-user tar -xf "${_rootfs}" -C "${_sandbox}" || true
        else
            tar -xf "${_rootfs}" --no-same-owner -C "${_sandbox}" || true
        fi
        rm -f "${_rootfs}"

        log "Checking the extracted rootfs"
        for _p in usr/local/freesurfer/*/bin/recon-all opt/ANTs/bin/antsRegistration \
                  opt/fsl/bin/fslmaths opt/mrtrix3/bin/mrcalc \
                  opt/KUL_VBG/KUL_VBG.sh opt/FastSurfer/run_fastsurfer.sh entrypoint.sh; do
            # shellcheck disable=SC2086  # glob on the FS version directory is intended
            compgen -G "${_sandbox}/${_p}" >/dev/null \
                || { rm -rf "${_work}"; die "flattened rootfs is incomplete: missing ${_p}"; }
        done

        # docker export carries only the filesystem — ENV/ENTRYPOINT/CMD are
        # image metadata and are lost, so they are restored here.
        cat > "${_work}/from_sandbox.def" <<EOF
Bootstrap: localimage
From: ${_sandbox}

%labels
    Maintainer   radwanphd@gmail.com
    Software     KUL_VBG
    Version      ${IMAGE_TAG}
    Source       https://github.com/KUL-Radneuron/KUL_VBG
    BuildRoute   docker export -> sandbox -> SIF (apptainer OCI unpack workaround)

%environment
    export KUL_VBG_FREESURFER_HOME=/usr/local/freesurfer/8.2.0
    export SUBJECTS_DIR=\${SUBJECTS_DIR:-/data}
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

%runscript
    exec /entrypoint.sh "\$@"
EOF

        log "Packing the sandbox into ${SIF_PATH}"
        "${APPTAINER_BIN}" build --force "${SIF_PATH}" "${_work}/from_sandbox.def" \
            || { rm -rf "${_work}"; die "apptainer build from sandbox failed"; }
        rm -rf "${_work}"
        ok "SIF built via the flatten route"
    fi

    ok "wrote ${SIF_PATH}  ($(du -h "${SIF_PATH}" | cut -f1))"

    # --cleanenv matters for the check itself, not just at run time: without it
    # apptainer merges the host environment and sources the user's ~/.profile
    # inside the container, which on a neuroimaging workstation refers to host
    # paths that do not exist here. That noise can mask a real problem.
    log "Verifying the SIF"
    # shellcheck disable=SC2016  # single quotes intended: this expands inside the container
    "${APPTAINER_BIN}" exec --cleanenv "${SIF_PATH}" bash -lc '
        set -e
        for t in recon-all mri_synthstrip mri_synthseg segment_subregions \
                 mri_segment_hypothalamic_subunits antsRegistration Atropos ImageMath \
                 antsRegistrationSyN.sh mrcalc mrstats labelconvert mrhistmatch \
                 fslmaths fslstats fslreorient2std imcp \
                 run_fastsurfer.sh KUL_VBG.sh KUL_VBG_multiparc.sh; do
            command -v "$t" >/dev/null || { echo "MISSING: $t"; exit 1; }
        done
        python3 -c "import numpy, scipy, nibabel, matplotlib" \
            || { echo "MISSING: python QC dependencies"; exit 1; }
        /opt/FastSurfer/.venv/bin/python -c "import torch" \
            || { echo "MISSING: FastSurfer torch runtime"; exit 1; }
        test -n "$(find /opt/FastSurfer/checkpoints -name "*.pkl" 2>/dev/null | head -1)" \
            || { echo "MISSING: FastSurfer checkpoints"; exit 1; }
        test -s /opt/KUL_VBG/atlasses/New/atlases/glasser/lh.HCPMMP1.annot \
            || { echo "MISSING: Glasser annot data"; exit 1; }
        # Must actually RUN, not merely exist. The unpatched FreeSurfer 8.2.0
        # .deb ships this binary but only the cpython-312 samseg bindings, so it
        # dies on import under the bundled Python 3.8. A command -v check passes
        # while every subregion segmentation (hippocampal, amygdala, thalamic,
        # brainstem) silently fails at run time. That is how a broken image
        # shipped once. NOTE: no apostrophes in here, this whole block is inside
        # a single-quoted bash -lc string.
        segment_subregions --help >/dev/null 2>&1 \
            || { echo "BROKEN: segment_subregions does not run (fs820 patch missing?)"; exit 1; }
    ' || die "SIF verification failed — the SIF is not usable"

    # The licence gate must REFUSE to run when nothing is bound. If this stops
    # working, users get an opaque FreeSurfer failure hours into a run instead
    # of an immediate, actionable error.
    if "${APPTAINER_BIN}" run --cleanenv "${SIF_PATH}" KUL_VBG.sh -h >/dev/null 2>&1; then
        die "SIF ran without a FreeSurfer licence — the licence gate is broken"
    fi
    ok "SIF verified (tools, torch, checkpoints, atlases, licence gate)"

    cat <<EOF

${C_HEAD}Next steps${C_OFF}
  Copy the SIF to the server:
      scp ${SIF_PATH} user@server:/path/to/containers/

  Run it there (no root needed):
      apptainer run --nv \\
          -B /path/to/license.txt:/licence/license.txt \\
          -B "\$PWD":/data \\
          /path/to/containers/$(basename "${SIF_PATH}") \\
          KUL_VBG.sh -S PAT001 -a /data/T1w.nii.gz -l /data/lesion.nii.gz \\
                     -z T1 -o /data/VBG_out -B 1 -P 1 -n 8 -v
EOF
fi

# ── Push ─────────────────────────────────────────────────────────────────────
if [[ ${DO_PUSH} -eq 1 ]]; then
    remote="${REGISTRY}/${IMAGE_NAME}:${IMAGE_TAG}"
    log "Pushing ${remote}"
    warn "This publishes the image to a public registry. Confirm you are logged in"
    warn "as the right account ('docker login') and that no licence file or patient"
    warn "data ever entered the build context."
    docker tag "${IMAGE_REF}" "${remote}"
    docker push "${remote}" || die "docker push failed"
    ok "pushed ${remote}"
fi

ok "done"
