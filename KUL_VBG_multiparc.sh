#!/bin/bash
# KUL_VBG_multiparc.sh
#
# Standalone multi-scale parcellation for any completed FreeSurfer recon-all.
# Runs the same atlas set as KUL_VBG.sh -M:
#   - Lausanne2018 scales 1-5  (mri_surf2surf + mri_aparc2aseg + MSBP ID remap)
#   - Glasser HCP-MMP1         (mri_surf2surf + mri_aparc2aseg)
#   - Thalamic nuclei          (segment_subregions thalamus,       FS 8+)
#   - Brainstem substructures  (segment_subregions brainstem,      FS 8+)
#   - Hippo/amygdala subregions(segment_subregions hippo-amygdala, FS 8+)
#   - Hypothalamic subunits    (mri_segment_hypothalamic_subunits, FS 7.2+)
#
# Requirements:
#   FreeSurfer 8+ on PATH (FREESURFER_HOME set), Python 3 + nibabel
#
# Usage:
#   KUL_VBG_multiparc.sh -S <subject_id> -f <fs_subjects_dir> [-n <threads>] [-v] [-O -l <lesion_mask>]
#
# Example:
#   KUL_VBG_multiparc.sh -S Patient01 -f /data/FS_output -n 16 -O -l /data/Patient01_lesion.nii.gz

version="0.3 — 2026-07-01"

# ── Defaults ──────────────────────────────────────────────────────────────────
ncpu=8
O_flag=0
verbose=0
subj=""
fs_dir="$(pwd)"
lesion_mask=""
t1_input=""

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
atlases_dir="${script_dir}/atlasses/New/atlases"
lausanne_dir="${atlases_dir}/lausanne2008"
glasser_dir="${atlases_dir}/glasser"
remap_py="${lausanne_dir}/remap_lausanne_to_msbp.py"

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
KUL_VBG_multiparc.sh v${version}

Usage: $(basename "$0") -S <subject_id> [options]

Required:
  -S <subject_id>  FreeSurfer subject ID (exact folder name under -f)

Options:
  -f <fs_dir>     FreeSurfer subjects directory (default: current directory)
  -T <T1.nii.gz>  T1-weighted image — runs recon-all if not already done
  -n <n>          Number of threads (default: ${ncpu})
  -v              Verbose: echo command output to the terminal as well as the log
  -O              Generate lesion overlap reports (requires -l)
  -l <file>       Binary lesion mask NIfTI — required with -O
  -h              Show this help

Outputs (written to <fs_subjects_dir>/sub-<subject_id>/mri/):
  lausanne2018.scale{1-5}+aseg.mgz
  HCPMMP1+aseg.mgz
  ThalamicNuclei.FSvoxelSpace.mgz        (FS 8+ only)
  brainstemSsLabels.FSvoxelSpace.mgz     (FS 8+ only)
  lh/rh.hippoAmygLabels.FSvoxelSpace.mgz (FS 8+ only)
  hypothalamic_subunits.v1.mgz           (FS 7.2+ only)

Done files (in <fs_subjects_dir>/sub-<subject_id>/scripts/):
  multiscale_parc.done       — Lausanne + Glasser complete
  thalamic_nuclei.done       — thalamus complete
  brainstem_subregions.done  — brainstem complete
  hippo_amygdala.done        — hippo-amygdala complete
  hypothalamic_subunits.done — hypothalamus complete
  Delete to force rerun of individual steps.

EOF
    exit 0
}

# ── Argument parsing ───────────────────────────────────────────────────────────
while getopts ":S:f:T:n:l:vOh" opt; do
    case $opt in
        S) subj="${OPTARG}" ;;
        f) fs_dir="$OPTARG" ;;
        T) t1_input="$OPTARG" ;;
        n) ncpu="$OPTARG" ;;
        l) lesion_mask="$OPTARG" ;;
        v) verbose=1 ;;
        O) O_flag=1 ;;
        h) usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
       \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
    esac
done

# ── Preflight checks ───────────────────────────────────────────────────────────
[[ -z "$subj" ]] && { echo "ERROR: -S <subject_id> required"; exit 1; }

subj_dir="${fs_dir}/${subj}"
if [[ ! -d "$subj_dir" ]] && [[ -z "${t1_input}" ]]; then
    echo "ERROR: subject dir not found: ${subj_dir}"
    echo "       Supply -T <T1.nii.gz> to run recon-all and create it."
    exit 1
fi

scripts_dir="${subj_dir}/scripts"

if [[ ! -f "${scripts_dir}/recon-all.done" ]] && [[ -z "${t1_input}" ]]; then
    echo "ERROR: recon-all.done not found in ${scripts_dir}"
    echo "       Either point -f/-S to a completed FreeSurfer directory, or"
    echo "       supply -T <T1.nii.gz> to run recon-all first."
    exit 1
fi

if [[ -n "${t1_input}" ]] && [[ ! -f "${t1_input}" ]]; then
    echo "ERROR: T1 image not found: ${t1_input}"; exit 1
fi

[[ -z "$FREESURFER_HOME" ]] && {
    echo "ERROR: FREESURFER_HOME not set — source FreeSurfer setup first"; exit 1
}

python3 -c "import nibabel, numpy" 2>/dev/null || {
    echo "ERROR: Python packages nibabel and numpy required (MSBP remap + overlap reports)"; exit 1
}
[[ ! -f "${remap_py}" ]] && { echo "ERROR: remap script not found: ${remap_py}"; exit 1; }

if [[ ${O_flag} -eq 1 ]]; then
    [[ -z "${lesion_mask}" ]] && { echo "ERROR: -l <lesion_mask> required with -O"; exit 1; }
    [[ ! -f "${lesion_mask}" ]] && { echo "ERROR: lesion mask not found: ${lesion_mask}"; exit 1; }
fi

_log_ts="$(date +%Y%m%d_%H%M%S)"
if [[ -d "${subj_dir}" ]]; then
    mkdir -p "${scripts_dir}"
    log_file="${scripts_dir}/KUL_VBG_multiparc_${_log_ts}.log"
else
    log_file="${fs_dir}/KUL_VBG_multiparc_${subj}_${_log_ts}.log"
fi

log() { echo "$*" | tee -a "${log_file}"; }

# ── Task runner ────────────────────────────────────────────────────────────────
# run: fatal on failure
run() {
    log "  "
    log "  $*"
    log "  Started @ $(date "+%Y-%m-%d_%H-%M-%S")"
    if [[ ${verbose} -eq 1 ]]; then
        eval "$*" 2>&1 | tee -a "${log_file}"
    else
        eval "$*" >> "${log_file}" 2>&1
    fi
    local _rc=${PIPESTATUS[0]}
    if [[ ${_rc} -eq 0 ]]; then
        log "  Success"
    else
        log "  Fail (exit ${_rc})"
        exit 1
    fi
    log "  Finished @ $(date "+%Y-%m-%d_%H-%M-%S")"
    log "  "
}

# run_soft: non-fatal, returns exit code
run_soft() {
    log "  "
    log "  $*"
    log "  Started @ $(date "+%Y-%m-%d_%H-%M-%S")"
    if [[ ${verbose} -eq 1 ]]; then
        eval "$*" 2>&1 | tee -a "${log_file}"
    else
        eval "$*" >> "${log_file}" 2>&1
    fi
    local _rc=${PIPESTATUS[0]}
    if [[ ${_rc} -eq 0 ]]; then
        log "  Success"
    else
        log "  Fail (non-fatal, exit ${_rc})"
    fi
    log "  Finished @ $(date "+%Y-%m-%d_%H-%M-%S")"
    log "  "
    return ${_rc}
}

# ── Overlap report helper ──────────────────────────────────────────────────────
# Converts MGZ → NIfTI (FS conformed space, no resampling) then runs Python overlap.
# Python script handles any remaining shape mismatch via nibabel resampling.
overlap_report() {
    # $1 parc  $2 report.txt  $3 name  $4 lut  $5 (optional) html path
    local _parc="$1" _report="$2" _name="$3" _lut="$4" _html="${5:-}"
    [[ ${O_flag} -eq 0 ]] && return 0
    [[ ! -f "${_parc}" ]] && { log "  [-O] ${_name} not found, skipping overlap"; return 0; }
    [[ ! -f "${_lut}" ]]  && { log "  [-O] LUT not found (${_lut}), skipping overlap"; return 0; }
    local _tmp _converted=0
    if [[ "${_parc}" == *.mgz ]]; then
        _tmp="${scripts_dir}/_ovl_tmp_$$.nii.gz"
        mri_convert "${_parc}" "${_tmp}" >>"${log_file}" 2>&1 \
            || { log "  [-O] mri_convert failed for ${_name}, skipping overlap"; return 0; }
        _converted=1
    else
        _tmp="${_parc}"
    fi
    local _html_arg="" _fallback_arg=""
    [[ -n "${_html}" ]] && _html_arg="--html ${_html} --subject ${subj}"
    [[ -f "${FREESURFER_HOME}/FreeSurferColorLUT.txt" ]] && \
        _fallback_arg="--lut-fallback ${FREESURFER_HOME}/FreeSurferColorLUT.txt"
    python3 "${script_dir}/KUL_lesion_overlap.py" \
        --parc "${_tmp}" --lesion "${lesion_mask}" \
        --out "${_report}" --name "${_name}" --lut "${_lut}" \
        ${_fallback_arg} ${_html_arg} \
        2>&1 | tee -a "${log_file}"
    log "  [-O] Overlap report → ${_report}"
    [[ ${_converted} -eq 1 ]] && rm -f "${_tmp}"
}

# ── Header ─────────────────────────────────────────────────────────────────────
log "========================================================"
log " KUL_VBG_multiparc.sh v${version}"
log " Subject:    ${subj}"
log " FS dir:     ${fs_dir}"
log " T1 input:   $( [[ -n "${t1_input}" ]] && echo "${t1_input}" || echo "(not provided — recon-all must already be done)" )"
log " Threads:    ${ncpu}"
log " Overlap:    $( [[ ${O_flag} -eq 1 ]] && echo "enabled (-O)" || echo disabled )"
log " Started:    $(date)"
log " Log:        ${log_file}"
log "========================================================"

export OMP_NUM_THREADS=${ncpu}
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${ncpu}

mri_dir="${subj_dir}/mri"

# ── FreeSurfer recon-all (optional, triggered by -T) ──────────────────────────
if [[ ! -f "${scripts_dir}/recon-all.done" ]]; then
    log ""
    log "recon-all.done absent — running FreeSurfer recon-all (this may take several hours)..."
    if [[ ! -f "${subj_dir}/mri/orig/001.mgz" ]]; then
        run "recon-all -s ${subj} -sd ${fs_dir} -i ${t1_input} -all -openmp ${ncpu}"
    else
        log "  Subject already initialized (orig/001.mgz exists) — continuing without -i"
        run "recon-all -s ${subj} -sd ${fs_dir} -all -openmp ${ncpu}"
    fi
    mkdir -p "${scripts_dir}"
    touch "${scripts_dir}/recon-all.done"
    if [[ "${log_file}" != "${scripts_dir}/"* ]]; then
        mv "${log_file}" "${scripts_dir}/KUL_VBG_multiparc_${_log_ts}.log"
        log_file="${scripts_dir}/KUL_VBG_multiparc_${_log_ts}.log"
    fi
    log "recon-all complete."
fi

# ── Lausanne2018 + Glasser ─────────────────────────────────────────────────────
if [[ -f "${scripts_dir}/multiscale_parc.done" ]]; then
    log "Lausanne/Glasser already done, skipping."
else
    # fsaverage symlink needed by mri_surf2surf
    [[ ! -d "${fs_dir}/fsaverage" ]] && \
        ln -sf "${FREESURFER_HOME}/subjects/fsaverage" "${fs_dir}/fsaverage"

    log ""
    log "Running Lausanne2018 parcellation (scales 1-5)..."

    for _scale in 1 2 3 4 5; do
        log "  Scale ${_scale}..."
        for _hemi in lh rh; do
            run "mri_surf2surf \
                --srcsubject fsaverage \
                --trgsubject ${subj} \
                --hemi ${_hemi} \
                --sval-annot ${lausanne_dir}/${_hemi}.lausanne2018.scale${_scale}.annot \
                --tval ${subj_dir}/label/${_hemi}.lausanne2018.scale${_scale}.annot \
                --sd ${fs_dir}"
        done
        _raw_parc="${mri_dir}/lausanne2018.scale${_scale}+aseg_raw.mgz"
        _lut_file="${lausanne_dir}/label-L2018_desc-scale${_scale}_atlas_FreeSurferColorLUT.txt"
        run "mri_aparc2aseg \
            --s ${subj} \
            --sd ${fs_dir} \
            --annot lausanne2018.scale${_scale} \
            --o ${_raw_parc}"
        if [[ -f "${_lut_file}" ]]; then
            run "python3 ${remap_py} \
                --input    ${_raw_parc} \
                --lh_annot ${lausanne_dir}/lh.lausanne2018.scale${_scale}.annot \
                --rh_annot ${lausanne_dir}/rh.lausanne2018.scale${_scale}.annot \
                --lut      ${_lut_file} \
                --output   ${mri_dir}/lausanne2018.scale${_scale}+aseg.mgz"
            rm -f "${_raw_parc}"
        else
            log "  WARNING: no MSBP LUT for scale ${_scale} — IDs not remapped"
            mv "${_raw_parc}" "${mri_dir}/lausanne2018.scale${_scale}+aseg.mgz" \
                || { log "ERROR: mv of ${_raw_parc} failed — NOT writing multiscale_parc.done"; exit 1; }
        fi
    done

    log ""
    log "Running Glasser HCP-MMP1 parcellation..."
    for _hemi in lh rh; do
        run "mri_surf2surf \
            --srcsubject fsaverage \
            --trgsubject ${subj} \
            --hemi ${_hemi} \
            --sval-annot ${glasser_dir}/${_hemi}.HCPMMP1.annot \
            --tval ${subj_dir}/label/${_hemi}.HCPMMP1.annot \
            --sd ${fs_dir}"
    done
    run "mri_aparc2aseg \
        --s ${subj} \
        --sd ${fs_dir} \
        --annot HCPMMP1 \
        --o ${mri_dir}/HCPMMP1+aseg.mgz"

    touch "${scripts_dir}/multiscale_parc.done"
    log "Lausanne/Glasser parcellation complete."
fi

# ── Overlap reports for Lausanne + Glasser ─────────────────────────────────────
if [[ ${O_flag} -eq 1 ]]; then
    _lesion_html="${scripts_dir}/${subj}_lesion_overlap_report.html"
    [[ -f "${_lesion_html}" ]] && rm -f "${_lesion_html}"
    for _scale in 1 2 3 4 5; do
        overlap_report \
            "${mri_dir}/lausanne2018.scale${_scale}+aseg.mgz" \
            "${scripts_dir}/lausanne_scale${_scale}_lesion_overlap.txt" \
            "Lausanne2018_scale${_scale}" \
            "${script_dir}/share/luts/lausanne_scale${_scale}_lut.txt" \
            "${_lesion_html}"
    done
    overlap_report \
        "${mri_dir}/HCPMMP1+aseg.mgz" \
        "${scripts_dir}/glasser_lesion_overlap.txt" \
        "Glasser_HCP-MMP1" \
        "${script_dir}/share/luts/glasser_lut.txt" \
        "${_lesion_html}"
fi

# ── Subcortical subsegmentations (FS 8+) ──────────────────────────────────────
# Each step has its own done file and is non-fatal (run_soft).
if ! command -v segment_subregions &>/dev/null; then
    log ""
    log "WARNING: segment_subregions not in PATH — FreeSurfer 8+ required for"
    log "         thalamic / brainstem / hippo-amygdala subsegmentations; skipping."
else
    # --- Thalamic nuclei ---
    if [[ ! -f "${scripts_dir}/thalamic_nuclei.done" ]]; then
        log ""
        log "Running thalamic nuclei segmentation (segment_subregions, FS 8+)..."
        if run_soft "segment_subregions thalamus --cross ${subj} --sd ${fs_dir} --threads ${ncpu}"; then
            touch "${scripts_dir}/thalamic_nuclei.done"
        else
            log "WARNING: thalamic segmentation failed — will retry on next run"
        fi
    else
        log "Thalamic nuclei already done, skipping."
    fi

    if [[ ${O_flag} -eq 1 ]]; then
        overlap_report \
            "${mri_dir}/ThalamicNuclei.FSvoxelSpace.mgz" \
            "${scripts_dir}/thalamic_nuclei_lesion_overlap.txt" \
            "ThalamicNuclei" \
            "${script_dir}/share/luts/thalamic_nuclei_lut.txt" \
            "${_lesion_html}"
    fi

    # --- Brainstem substructures ---
    if [[ ! -f "${scripts_dir}/brainstem_subregions.done" ]]; then
        log ""
        log "Running brainstem substructure segmentation (segment_subregions, FS 8+)..."
        if run_soft "segment_subregions brainstem --cross ${subj} --sd ${fs_dir} --threads ${ncpu}"; then
            touch "${scripts_dir}/brainstem_subregions.done"
        else
            log "WARNING: brainstem segmentation failed — will retry on next run"
        fi
    else
        log "Brainstem substructures already done, skipping."
    fi

    if [[ ${O_flag} -eq 1 ]]; then
        overlap_report \
            "${mri_dir}/brainstemSsLabels.FSvoxelSpace.mgz" \
            "${scripts_dir}/brainstem_lesion_overlap.txt" \
            "BrainstemSubstructures" \
            "${script_dir}/share/luts/brainstem_lut.txt" \
            "${_lesion_html}"
    fi

    # --- Hippocampal/amygdala subregions ---
    if [[ ! -f "${scripts_dir}/hippo_amygdala.done" ]]; then
        log ""
        log "Running hippocampal/amygdala subregion segmentation (segment_subregions, FS 8+)..."
        if run_soft "segment_subregions hippo-amygdala --cross ${subj} --sd ${fs_dir} --threads ${ncpu}"; then
            touch "${scripts_dir}/hippo_amygdala.done"
        else
            log "WARNING: hippo-amygdala segmentation failed — will retry on next run"
        fi
    else
        log "Hippo-amygdala segmentation already done, skipping."
    fi

    if [[ ${O_flag} -eq 1 ]]; then
        overlap_report \
            "${mri_dir}/lh.hippoAmygLabels.FSvoxelSpace.mgz" \
            "${scripts_dir}/lh_hippo_amyg_lesion_overlap.txt" \
            "HippoAmyg-LH" \
            "${script_dir}/share/luts/hippo_amygdala_lut.txt" \
            "${_lesion_html}"
        overlap_report \
            "${mri_dir}/rh.hippoAmygLabels.FSvoxelSpace.mgz" \
            "${scripts_dir}/rh_hippo_amyg_lesion_overlap.txt" \
            "HippoAmyg-RH" \
            "${script_dir}/share/luts/hippo_amygdala_lut.txt" \
            "${_lesion_html}"
    fi
fi

# ── Hypothalamic subunits (FS 7.2+) ───────────────────────────────────────────
if ! command -v mri_segment_hypothalamic_subunits &>/dev/null; then
    log ""
    log "WARNING: mri_segment_hypothalamic_subunits not in PATH — skipping hypothalamus."
else
    if [[ ! -f "${scripts_dir}/hypothalamic_subunits.done" ]]; then
        log ""
        log "Running hypothalamic subunit segmentation..."
        if run_soft "mri_segment_hypothalamic_subunits --s ${subj} --sd ${fs_dir} --threads ${ncpu}"; then
            touch "${scripts_dir}/hypothalamic_subunits.done"
        else
            log "WARNING: hypothalamic segmentation failed — will retry on next run"
        fi
    else
        log "Hypothalamic subunits already done, skipping."
    fi

    if [[ ${O_flag} -eq 1 ]]; then
        overlap_report \
            "${mri_dir}/hypothalamic_subunits.v1.mgz" \
            "${scripts_dir}/hypothalamic_subunits_lesion_overlap.txt" \
            "HypothalamicSubunits" \
            "${script_dir}/share/luts/hypothalamic_subunits_lut.txt" \
            "${_lesion_html}"
        log "  [-O] HTML report → ${_lesion_html}"
    fi
fi

# ── Done ──────────────────────────────────────────────────────────────────────
log ""
log "========================================================"
log " KUL_VBG_multiparc.sh complete."
log " Outputs in ${mri_dir}/"
log " Finished: $(date)"
log "========================================================"
