#!/bin/bash
###############################################################################
# KUL_VBG container entrypoint — Docker and Apptainer.
#
# Deliberately written to work identically under both runtimes, which means it
# must never write anything inside the image: an Apptainer SIF is read-only.
# (The previous version copied the mounted FreeSurfer licence into
# $FREESURFER_HOME and into /opt/FastSurfer — that fails outright on a SIF.)
#
# Responsibilities, in order:
#   1. locate and validate the user-supplied FreeSurfer licence
#   2. export the environment KUL_VBG.sh expects
#   3. check /tmp is writable (the FS 8.2.0 mris_register workaround needs it)
#   4. report GPU visibility
#   5. exec whatever the user asked for
###############################################################################
set -euo pipefail

_die() { echo "ERROR: $*" >&2; exit 1; }

###############################################################################
# 1. FreeSurfer licence
#
# Not baked into the image on purpose — FreeSurfer's licence is per-user and
# redistributing it inside a shared image would violate its terms. Get one free
# at https://surfer.nmr.mgh.harvard.edu/registration.html
#
# $FREESURFER_HOME/license.txt is a symlink to /licence/license.txt (see the
# Dockerfile), so the canonical thing to bind is /licence/license.txt. A few
# other common paths are accepted too, as is FS_LICENSE if the user set it.
###############################################################################
_lic=""
for _c in \
    "${FS_LICENSE:-}" \
    /licence/license.txt \
    /license/license.txt \
    "${FREESURFER_HOME:-}/license.txt" \
    "${FREESURFER_HOME:-}/.license" \
    /opt/freesurfer/license.txt
do
    # -s not -f: the image ships an EMPTY placeholder at /licence/license.txt so
    # that the bind target exists on clusters without overlay support. An empty
    # file means "nothing was bound here", not "here is a licence" — accepting
    # it would surface as a confusing FreeSurfer failure much later in the run.
    if [[ -n "${_c}" && -s "${_c}" ]]; then
        _lic="${_c}"
        break
    fi
done

if [[ -z "${_lic}" ]]; then
    cat >&2 <<'EOF'
ERROR: no FreeSurfer licence found.

KUL_VBG needs a FreeSurfer licence at run time. It is deliberately NOT included
in this image — the licence is per-user and cannot be redistributed. Register
for free at: https://surfer.nmr.mgh.harvard.edu/registration.html

Then bind it in at /licence/license.txt:

  Docker:
    docker run --rm -it \
      -v /path/to/license.txt:/licence/license.txt:ro \
      -v "$PWD":/data \
      kul_vbg:2.0 KUL_VBG.sh -h

  Apptainer:
    apptainer run \
      -B /path/to/license.txt:/licence/license.txt \
      -B "$PWD":/data \
      KUL_VBG_2.0.sif KUL_VBG.sh -h
EOF
    exit 2
fi

export FS_LICENSE="${_lic}"

###############################################################################
# 2. Environment
#
# The Dockerfile sets all of this via ENV, so it is already correct for a plain
# `docker run`. Re-asserting it here matters for Apptainer: unlike Docker,
# Apptainer by default merges the HOST environment into the container, so a
# user whose login shell exports its own FSLDIR / ANTSPATH / FREESURFER_HOME
# (very common on a neuroimaging server — which is exactly where this image is
# meant to run) would otherwise have the host's paths silently shadow the
# container's, producing "command not found" or, worse, a mix of two FreeSurfer
# installs. Pinning them here makes the container authoritative.
###############################################################################
export FREESURFER_HOME="${KUL_VBG_FREESURFER_HOME:-/usr/local/freesurfer/8.2.0}"
export FSLDIR=/opt/fsl
export FSLOUTPUTTYPE=NIFTI_GZ
export ANTSPATH=/opt/ANTs/bin
export FASTSURFER_HOME=/opt/FastSurfer
export PATH="/opt/KUL_VBG:/opt/vbg-pyenv/bin:${ANTSPATH}:/opt/mrtrix3/bin:${FSLDIR}/bin:${FSLDIR}/share/fsl/bin:${FREESURFER_HOME}/bin:${FASTSURFER_HOME}:/usr/local/bin:/usr/bin:/bin"

# SetUpFreeSurfer.sh sets a handful of variables beyond FREESURFER_HOME
# (FSFAST_HOME, MNI_DIR, MINC paths...). It is noisy, references unset
# variables, and on some versions returns non-zero harmlessly — hence the
# `set +u` guard and the discarded output/status. Everything that actually
# matters is asserted explicitly above and below.
if [[ -f "${FREESURFER_HOME}/SetUpFreeSurfer.sh" ]]; then
    set +u
    # shellcheck disable=SC1091  # path only exists inside the image
    source "${FREESURFER_HOME}/SetUpFreeSurfer.sh" >/dev/null 2>&1 || true
    set -u
fi

# SUBJECTS_DIR defaults to the bind point rather than anywhere inside the
# image, so recon-all output lands on the user's filesystem.
export SUBJECTS_DIR="${SUBJECTS_DIR:-/data}"

###############################################################################
# 3. Writable /tmp
#
# KUL_VBG.sh works around a confirmed FreeSurfer 8.2.0 mris_register buffer
# overflow on long SUBJECTS_DIR paths by redirecting SUBJECTS_DIR through a
# short symlink under /tmp (KUL_shorten_SUBJECTS_DIR). Docker always gives a
# writable /tmp; Apptainer normally does too, but some sites bind it read-only
# or to a quota'd location. Failing here with a clear message beats crashing
# hours into recon-all.
###############################################################################
if ! : > /tmp/.kul_vbg_write_test 2>/dev/null; then
    _die "/tmp is not writable inside the container.
KUL_VBG needs it for the FreeSurfer 8.2.0 mris_register workaround
(SUBJECTS_DIR is redirected through a short symlink under /tmp).
Bind a writable directory over it, e.g.  -B /scratch/\$USER/tmp:/tmp"
fi
rm -f /tmp/.kul_vbg_write_test

###############################################################################
# 4. GPU visibility — informational only
#
# VBG runs fine CPU-only; FastSurfer (-P 2/3) is just much slower. Reporting
# this up front avoids the common "why has this been running for six hours"
# support question, and catches the frequent Apptainer mistake of forgetting
# --nv (without it the host driver is not injected and torch sees no GPU even
# on a GPU node).
###############################################################################
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
    echo "INFO: GPU visible to the container — $(nvidia-smi -L 2>/dev/null | head -1)"
else
    echo "INFO: no GPU visible; FastSurfer (-P 2/3) will run on CPU (slower but correct)."
    echo "      If this machine HAS a GPU, you likely forgot 'docker run --gpus all' or 'apptainer --nv'."
fi

###############################################################################
# 5. Hand off
###############################################################################
if [[ $# -eq 0 ]]; then
    exec KUL_VBG.sh -h
fi

# Bare `bash`/`sh` gets an interactive login shell so /etc/profile.d is picked
# up; anything else runs directly with the environment set above.
case "$1" in
    bash|sh|/bin/bash|/bin/sh)
        shift
        exec /bin/bash --login "$@"
        ;;
esac

exec "$@"
