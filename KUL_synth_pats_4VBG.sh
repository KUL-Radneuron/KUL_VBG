#!/bin/bash

set -x

# @ Ahmed Radwan -> ahmed.radwan@kuleuven.be, radwanphd@gmail.com
# @ Stefan Sunaert -> stefan.sunaert@kuleuven.be
# KUL_synth_cohort_gen.sh

# v=0.25 - 11092021

cwd="$(pwd)"

# function Usage
function Usage {

cat <<USAGE

    `basename $0` Runs whole brain TCK segmentation using an input config file

    Usage:

    `basename $0` -P /path_to/patients_dir -H /path_to/healthy_controls_dir -M /path_to/healthy_controls_dir -c /path_to/get_lesions_config_file.txt -o /path_to/Output_dir -n 6 -R 1

    Examples:

    `basename $0` -P /path_to/patients_dir -H /path_to/healthy_controls_dir -M /path_to/healthy_controls_dir -c /path_to/get_lesions_config_file.txt -o /path_to/Output_dir -n 6 -R 1

    Purpose:

    This workflow generates synthetic patients images, matching mass effect in a patient population without the lesion + synthetic lesioned brains
    We use a group of patients with a focal lesion, the lesion mask, and VBG lesion filled native images, in combination with healthy volunteer whole head T1 images
    The generated dataset consists of the product of patients X healthy volunteers X 2, as every combination is generated once only with deformation to match mass effect
    and once with the lesion. Patient and Healthy volunteer T1-weighted images should be acquired from similar scanners to ensure similar signal patterns. 
    An alternative to this is to use rescaled T1-weighted images e.g. the nu.mgz images of FreeSurfer's recon-all after conversion to .nii.gz

    Required arguments:

    -P:  Full path to patients' directory (this should contain a pair of native space VBG filled and original lesioned T1 whole head images, in 1 subfolder per patient)
    -H:  Full path to healthy controls' directory (this should contain all healthy volunteer images, each in a separate subfolder)
    -M:  Full path to patients' lesion masks (this should contain all patients lesion masks in T1 space, each in a separate subfolder)
    -c:  Full path and file name of config file containing all patients in column 1 and all controls in column 2, columns should be separated by commas
    * subfolders should be named after the subjects, e.g. sub-PT001, sub-PT002, etc.

    Optional arguments:

    -R:  Specify registration approach used (1: custom ants registration as used for the VBG paper, 2: antsRegSyN.sh, 3: antsRegSyNQuick.sh)
    -o:  Full path to output directory
    -n:  Number of cpu for parallelisation (default is 6)
    -h:  Prints help menu

    Example config file:

    PT001,HV01
    PT002,HV02
    PT003,HV03


USAGE

    exit 1
}

# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
# this works for ANTsX scripts and FS
# Set required options
P_flag=0
H_flag=0
M_flag=0
o_flag=0
c_flag=0
n_flag=0
R_flag=0

if [ "$#" -lt 5 ]; then
    Usage >&2
    exit 1

else

    while getopts "P:H:M:R:o:c:n:h" OPT; do

        case $OPT in
        P) #Patients
            P_flag=1
            pats_d=$OPTARG
        ;;
        H) #Controls
            H_flag=1
            cons_d=$OPTARG
        ;;
        M) #Masks
            M_flag=1
            masks_d=$OPTARG
        ;;
        R) #Registration approach
            R_flag=1
            reg=$OPTARG
        ;;
        c) #config file
            c_flag=1
            conf_f=$OPTARG
        ;;
        o) #output
            o_flag=1
            out_dir=$OPTARG
        ;;
        n) #parallel
            n_flag=1
            ncpu=$OPTARG
        ;;
        h) #help
            Usage >&2
            exit 0
        ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo
            Usage >&2
            exit 1
        ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            echo
            Usage >&2
            exit 1
        ;;
        esac

    done

fi

# also need to make sure we find the pats and HVs in the config file
# config file
srch_conf_str=($(basename ${conf_f})) ; conf_dir=($(dirname ${conf_f}))
srch_conf_c=($(find ${conf_dir} -type f | grep  ${srch_conf_str}))

# Pats dir
srch_Pdir_str=($(basename ${pats_d})) ; Pats_dir=($(dirname ${pats_d}))
srch_Pdir_c=($(find ${Pats_dir} -type d | grep  ${srch_Pdir_str}))

# Healthy controls dir
srch_Hdir_str=($(basename ${cons_d})) ; HC_dir=($(dirname ${cons_d}))
srch_Hdir_c=($(find ${HC_dir} -type d | grep  ${srch_Hdir_str}))

# Pats lesion masks dir
srch_Ldir_str=($(basename ${masks_d})) ; LM_dir=($(dirname ${masks_d}))
srch_Ldir_c=($(find ${LM_dir} -type d | grep  ${srch_Ldir_str}))

if [[ ${P_flag} -eq 0 ]] || [[ ${H_flag} -eq 0 ]] || [[ ${M_flag} -eq 0 ]] || [[ ${o_flag} -eq 0 ]] || [[ ${c_flag} -eq 0 ]]; then

    echo "incorrect input arguments, quitting"
    exit 2

else

    if [[ -z "${srch_Pdir_c}" ]]; then

        echo
        echo " Incorrect path to Patients' brain images dir, please check the path and name "
        echo
        exit 2

    fi

    if [[ -z "${srch_Hdir_c}" ]]; then

        echo
        echo " Incorrect path to the Healthy control brain images dir, please check the path and name "
        echo
        exit 2

    fi

    if [[ -z "${srch_Ldir_c}" ]]; then

        echo
        echo " Incorrect path to the lesion masks dir, please check the path and name "
        echo
        exit 2

    fi

    if [[ -z "${srch_conf_c}" ]]; then

        echo
        echo " Incorrect config file, please check the path and name "
        echo
        exit 2

    fi

    echo "inputs are -P ${pats_d} -H ${cons_d} -M ${masks_d}"

fi

# exit 2
# deal with nthreads output dir, prep dir and log file

if [[ "$n_flag" -eq 0 ]]; then

	ncpu=6

	echo " -n flag not set, using default 8 threads. "

else

	echo " -n flag set, using " ${ncpu} " threads."

fi

FSLPARALLEL=$ncpu; export FSLPARALLEL
OMP_NUM_THREADS=$ncpu; export OMP_NUM_THREADS

#

if [[ "$R_flag" -eq 0 ]]; then

	reg=1

	echo " -R flag not set, using custom ANTs registration. "
    echo " you have chosen to use customized ANTs registration, this takes a while" | tee -a ${prep_log}

elif [[ "$R_flag" -eq 1 ]]; then

	echo " -R flag set, using custom ANTs registration "
    echo " you have chosen to use customized ANTs registration, this takes a while" | tee -a ${prep_log}

elif [[ "$R_flag" -eq 2 ]]; then

    echo " -R flag set, using custom antsRegistrationSyN.sh "
    echo " you have chosen to use antsRegistrationSyN.sh " | tee -a ${prep_log}

elif [[ "$R_flag" -eq 3 ]]; then

    echo " -R flag set, using custom antsRegistrationSyNQuick.sh "
    echo " you have chosen to use antsRegistrationSyNQuick.sh, this is rather quick and dirty " | tee -a ${prep_log}

fi

# timestamp
start=$(date +%s)
d=$(date "+%Y-%m-%d_%H-%M-%S")

# handle the dirs
cwd=$(pwd)

cd ${cwd}

# handle output and processing dirs

if [[ "$o_flag" -eq 1 ]]; then

    output_d="${out_dir}"

else

    output_d="${cwd}/KUL_VBG_synth_pats_output"

fi

# output sub-dirs

out_SME="${output_d}/Synthetic_mass_effect_patients"

out_SP="${output_d}/Synthetic_lesioned_patients"

int="${output_d}/temp_d"

# make your dirs

mkdir -p ${output_d} >/dev/null 2>&1

mkdir -p ${out_SME} >/dev/null 2>&1

mkdir -p ${out_SP} >/dev/null 2>&1

mkdir -p ${int} >/dev/null 2>&1

# make your log file

prep_log="${output_d}/KUL_synth_pats_${d}.txt";

if [[ ! -f ${prep_log} ]] ; then

    touch ${prep_log}

else

    echo "${prep_log} already created"

fi

# set mrtrix tmp dir to tmpo_d

rm -rf ${int}/tmp_ims_*

tmpo_d=($(find ${int} -type d -name *"tmp_ims_"*))

if [[ -z ${tmpo_d} ]]; then

    tmpo_d="${int}/tmp_ims_${d}"

fi

mkdir -p "${tmpo_d}" >/dev/null 2>&1

export MRTRIX_TMPFILE_DIR="${tmpo_d}"

# report pid

processId=$(ps -ef | grep 'ABCD' | grep -v 'grep' | awk '{ printf $2 }')
echo $processId

# start by reading the config file

# tck_lst1=($(cat  ${conf_f}))

IFS=$'\n' read -d '' -r -a all_subs < ${conf_f}

for i in ${!all_subs[@]}; do

    if [[ ${all_subs[$i]} == *"#"* ]]; then

        PTs[$i]="none"
        HVs[$i]="none"

    else

        # tck_list[$i]=${tck_lst1[$i]}
        PTs[$i]=$(echo ${all_subs[$i]} | cut -d ',' -f1)
        # test to make sure tck_list[$i] contains a string
        if [[ ${PTs[$i]} =~ ^[+-]?[0-9]+$ ]]; then 
            echo " there is a problem with config file, first column does not contain a string" 
            exit 2
        fi

        HVs[$i]=$(echo ${all_subs[$i]} | cut -d ',' -f2)
        if [[ ${HVs[$i]} =~ ^[+-]?[0-9]+$ ]]; then 
            echo " there is a problem with config file, second column does not contain a string" 
            exit 2
        fi

    fi

done

# tell the use what we found

echo "You have specified the following subjects, Patients: ${Pats[@]},  Healthy volunteers: ${HVs[@]}"

# define task_exec function

function task_exec {

    echo "-------------------------------------------------------------" | tee -a ${prep_log}

    echo ${task_in} | tee -a ${prep_log}

    echo " Started @ $(date "+%Y-%m-%d_%H-%M-%S")" | tee -a ${prep_log}

    eval ${task_in} 2>&1 | tee -a ${prep_log} &

    # echo " pid = $! basicPID = $BASHPID " | tee -a ${prep_log}

    echo " pid = $! " | tee -a ${prep_log}

    wait ${pid}

    sleep 5

    echo "exit status $?" | tee -a ${prep_log}

    # if [ $? -eq 0 ]; then

    #     echo Success >> ${prep_log}

    # else

    #     echo Fail >> ${prep_log}

    #     exit 1

    # fi

    echo " Finished @ $(date "+%Y-%m-%d_%H-%M-%S")" | tee -a ${prep_log}

    echo "-------------------------------------------------------------" | tee -a ${prep_log}

    echo "" | tee -a ${prep_log}

    unset task_in

}

# please cite us

echo " You are using the synthetic cohort generation workflow, part of the VBG package, please cite the following paper in your work: \
Virtual brain grafting: Enabling whole brain parcellation in the presence of large lesions \
Ahmed M. Radwan, Louise Emsell, Jeroen Blommaert, Andrey Zhylka, Silvia Kovacs, Tom Theys, Nico Sollmann, Patrick Dupont, Stefan Sunaert \
medRxiv 2020.09.30.20204701; doi: https://doi.org/10.1016/j.neuroimage.2021.117731" | tee -a ${prep_log}


# start of script

declare -a PT_d_in
declare -a PT_intd
declare -a HV_intd
declare -a PT_masks
declare -a PTinHV_int
declare -a PT_outad
declare -a PT_outbd

# pow=0

# pew=${ncpu}

nvd_cu=$(nvcc --version)

for pt in ${!PTs[@]}; do

    if [[ ! ${PTs[$pt]} == "none" ]]; then

        PT_masks[$pt]="${masks_d}/sub-${PTs[$pt]}/sub-${PTs[$pt]}_Lesion_mask.nii.gz"

        PT_intd[$pt]="${int}/sub-${PTs[$pt]}_proc"

        mkdir -p ${PT_intd[$pt]}

        # define T1s, BM based on the dir
        PT_d_in[$pt]="${pats_d}/sub-${PTs[$pt]}"

    fi

done

for im in ${!HVs[@]}; do

    if [[ ! ${HVs[$im]} == "none" ]]; then

        HV_intd[$im]="${int}/sub-${HVs[$im]}_proc"

        mkdir -p ${HV_intd[$im]}

    fi

done

for ia in ${!HVs[@]}; do

    for pti in ${!PTs[@]}; do

        if [[ ! ${PTs[$pti]} == "none" ]]; then

            PTinHV_int[$pti]="${HV_intd[$ia]}/sub-${PTs[$pti]}_in_${HVs[$ia]}"

            mkdir -p ${PTinHV_int[$pti]}

            PT_outad[$pti]="${out_SME}/sub-${PTs[$pti]}_in_${HVs[$ia]}"

            PT_outbd[$pti]="${out_SP}/sub-${PTs[$pti]}_in_${HVs[$ia]}"

            mkdir -p ${PT_outad[$pti]}

            mkdir -p ${PT_outbd[$pti]}

        fi
    
    done

done


for ptx in ${!PTs[@]}; do

    if [[ ! ${PTs[$ptx]} == "none" ]]; then

        # some modifications are needed here
        # let's extract the brain again for the patients
        # we use this brain and brain mask for the registration
        # This mask should be applied to the non-filled T1s also 

        if [[ ! -f "${PT_intd[$ptx]}/sub-${PTs[$ptx]}_T1w_filled_brain.nii.gz" ]]; then

            task_in="recon-all -i ${PT_d_in[$ptx]}/sub-${PTs[$ptx]}_T1w.nii.gz -s ${PTs[$ptx]}_temp \
            -sd ${PT_intd[$ptx]} -openmp ${ncpu} -parallel -autorecon1"

            task_exec

            # convert to nii from the initial autorecon1 output
            task_in="mri_convert -rl ${PT_d_in[$ptx]}/sub-${PTs[$ptx]}_T1w.nii.gz \
            ${PT_intd[$ptx]}/${PTs[$ptx]}_temp/mri/orig_nu.mgz ${PT_intd[$ptx]}/sub-${PTs[$ptx]}_T1_nu.nii.gz"

            task_exec

            if [[ -z ${nvd_cu} ]]; then
                
                task_in="hd-bet -i ${PT_d_in[$ptx]}/sub-${PTs[$ptx]}_T1_nat_filled.nii.gz \
                -o ${PT_intd[$ptx]}/sub-${PTs[$ptx]}_T1w_filled_brain -tta 0 -mode accurate -s 1 -device cpu"

                task_exec

            else

                task_in="hd-bet -i ${PT_d_in[$ptx]}/sub-${PTs[$ptx]}_T1_nat_filled.nii.gz \
                -o ${PT_intd[$ptx]}/sub-${PTs[$ptx]}_T1w_filled_brain -mode accurate -s 1"

                task_exec

            fi

        else

            echo " ${PT_intd[$ptx]}/sub-${PTs[$ptx]}_T1w_filled_brain already done"

        fi
          
        task_in="fslmaths ${PT_intd[$ptx]}/sub-${PTs[$ptx]}_T1_nu.nii.gz \
        -mas ${PT_intd[$ptx]}/sub-${PTs[$ptx]}_T1w_filled_brain_mask.nii.gz ${PT_intd[$ptx]}/sub-${PTs[$ptx]}_rT1w_brain.nii.gz"

        task_exec

        if [[ ! -f "${PT_intd[$ptx]}/sub-${PTs[$ptx]}_LM_dst_inv_BM.nii.gz" ]]; then 

            echo "Now working on making inv Lesion masks" | tee -a ${prep_log}

            task_in="fslmaths ${PT_masks[$ptx]} -dilM -s 2 -thr 0.3 -save \
            ${PT_intd[$ptx]}/sub-${PTs[$ptx]}_LM_dst.nii.gz -binv -mul \
            ${PT_intd[$ptx]}/sub-${PTs[$ptx]}_T1w_filled_brain_mask.nii.gz \
            ${PT_intd[$ptx]}/sub-${PTs[$ptx]}_LM_dst_inv_BM.nii.gz"

            task_exec

        else

            echo " ${PT_intd[$ptx]}/sub-${PTs[$ptx]} lesion mask derivatives already done" | tee -a ${prep_log}

        fi

    fi

done


# Loop over HVs

for i in ${!HVs[@]}; do

    if [[ ! ${HVs[$i]} == "none" ]]; then

        # extract HVs brains

        echo "Now working on HVs BETs" | tee -a ${prep_log}

        # look for CUDA installation
        # if not found use hd-bet in CPU mode
        # if found use hd-bet with cuda

        if [[ ! -f "${HV_intd[$i]}/sub-${HVs[$i]}_T1w_brain.nii.gz" ]]; then

            task_in="fslreorient2std ${cons_d}/sub-${HVs[$i]}/sub-${HVs[$i]}_T1w.nii.gz ${HV_intd[$i]}/sub-${HVs[$i]}_T1w_reori.nii.gz"

            task_exec

            task_in="recon-all -i ${HV_intd[$i]}/sub-${HVs[$i]}_T1w_reori.nii.gz -s ${HVs[$i]}_temp \
            -sd ${HV_intd[$i]} -openmp ${ncpu} -parallel -autorecon1"

            task_exec

            # convert to nii from the initial autorecon1 output
            task_in="mri_convert -rl ${HV_intd[$i]}/sub-${HVs[$i]}_T1w_reori.nii.gz \
            ${HV_intd[$i]}/${HVs[$i]}_temp/mri/orig_nu.mgz ${HV_intd[$i]}/sub-${HVs[$i]}_T1_nu.nii.gz"

            task_exec

            if [[ -z ${nvd_cu} ]]; then

                task_in="hd-bet -i ${HV_intd[$i]}/sub-${HVs[$i]}_T1w_reori.nii.gz \
                -o ${HV_intd[$i]}/sub-${HVs[$i]}_T1w_brain -tta 0 -mode accurate -s 1 -device cpu"

                task_exec
            else

                task_in="hd-bet -i ${HV_intd[$i]}/sub-${HVs[$i]}_T1w_reori.nii.gz \
                -o ${HV_intd[$i]}/sub-${HVs[$i]}_T1w_brain -mode accurate -s 1"

                task_exec

            fi


        else

            echo " ${HV_intd[$i]}/sub-${HVs[$i]}_T1_brain already done"

        fi

        task_in="fslmaths ${HV_intd[$i]}/sub-${HVs[$i]}_T1_nu.nii.gz \
        -mas ${HV_intd[$i]}/sub-${HVs[$i]}_T1w_brain_mask.nii.gz ${HV_intd[$i]}/sub-${HVs[$i]}_rT1w_brain.nii.gz"

        task_exec

        # loop over all pats per HV

        for p in ${!PTs[@]}; do

            if [[ ! ${PTs[$p]} == "none" ]]; then

                # # define T1s, BM based on the dir
                # PT_d_in[$p]="${pats_d}/sub-${PTs[$p]}"

                PTinHV_int[$p]="${HV_intd[$i]}/sub-${PTs[$p]}_in_${HVs[$i]}"

                # mkdir -p ${PTinHV_int[$p]}

                PT_outad[$p]="${out_SME}/sub-${PTs[$p]}_in_${HVs[$i]}"

                PT_outbd[$p]="${out_SP}/sub-${PTs[$p]}_in_${HVs[$i]}"

                # mkdir -p ${PT_outad[$p]}

                # mkdir -p ${PT_outbd[$p]}

                # some modifications are needed here
                # let's extract the brain again for the patients
                # we use this brain and brain mask for the registration
                # This mask should be applied to the non-filled T1s also 

                # exit 2

                if [[ ! -f "${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_Warped.nii.gz" ]]; then

                    # ants warping of filled pat T1 to HV T1

                    if [[ ${reg} -eq 1 ]]; then

                        task_in="export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${ncpu} ; antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 -u 1 \
                        --output [ ${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_,${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_Warped.nii.gz,${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_InverseWarped.nii.gz ] \
                        --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [ 0.005,0.995 ] \
                        -x [ ${HV_intd[$i]}/sub-${HVs[$i]}_T1w_brain_mask.nii.gz,${PT_intd[$p]}/sub-${PTs[$p]}_T1w_filled_brain_mask.nii.gz, NULL ] \
                        --initial-moving-transform [ ${HV_intd[$i]}/sub-${HVs[$i]}_rT1w_brain.nii.gz,${PT_intd[$p]}/sub-${PTs[$p]}_T1w_filled_brain.nii.gz,1 ] --transform Rigid[ 0.1 ] \
                        --metric MI[ ${HV_intd[$i]}/sub-${HVs[$i]}_rT1w_brain.nii.gz,${PT_intd[$p]}/sub-${PTs[$p]}_T1w_filled_brain.nii.gz,1,32,Regular,0.25 ] \
                        --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[ 0.1 ] \
                        --metric MI[ ${HV_intd[$i]}/sub-${HVs[$i]}_rT1w_brain.nii.gz,${PT_intd[$p]}/sub-${PTs[$p]}_T1w_filled_brain.nii.gz,1,64,Regular,0.5 ] \
                        --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform SyN[ 0.1,3,0 ] \
                        --metric CC[ ${HV_intd[$i]}/sub-${HVs[$i]}_rT1w_brain.nii.gz,${PT_intd[$p]}/sub-${PTs[$p]}_T1w_filled_brain.nii.gz,1,4 ] \
                        --convergence [ 200x100x75x25,1e-8,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --verbose 1"

                        task_exec

                    elif [[ ${reg} -eq 2 ]]; then

                        task_in="antsRegistrationSyN.sh -d 3 -f ${HV_intd[$i]}/sub-${HVs[$i]}_rT1w_brain.nii.gz -m ${PT_intd[$p]}/sub-${PTs[$p]}_T1w_filled_brain.nii.gz \
                        -x ${HV_intd[$i]}/sub-${HVs[$i]}_T1w_brain_mask.nii.gz,${PT_intd[$p]}/sub-${PTs[$p]}_T1w_filled_brain_mask.nii.gz -t s -n ${ncpu} -o ${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_"

                        task_exec

                    elif [[ ${reg} -eq 3 ]]; then

                        task_in="antsRegistrationSyNQuick.sh -d 3 -f ${HV_intd[$i]}/sub-${HVs[$i]}_rT1w_brain.nii.gz -m ${PT_intd[$p]}/sub-${PTs[$p]}_T1w_filled_brain.nii.gz \
                        -x ${HV_intd[$i]}/sub-${HVs[$i]}_T1w_brain_mask.nii.gz,${PT_intd[$p]}/sub-${PTs[$p]}_T1w_filled_brain_mask.nii.gz -t s -n ${ncpu} -o ${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_"

                        task_exec

                    fi

                fi

                # exit 2

                if [[ ! -f "${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_intmatched.nii.gz" ]]; then

                    task_in="antsApplyTransforms -d 3 -i ${HV_intd[$i]}/sub-${HVs[$i]}_T1_nu.nii.gz -o ${PTinHV_int[$p]}/sub-${HVs[$i]}_T1w_in_${PTs[$p]}_Warped1.nii.gz -r ${PT_d_in[$p]}/sub-${PTs[$p]}_T1_nat_filled.nii.gz \
                    -t [${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_0GenericAffine.mat,1] -t ${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_1InverseWarp.nii.gz -n LanczosWindowedSinc \
                    && antsApplyTransforms -d 3 -i ${HV_intd[$i]}/sub-${HVs[$i]}_T1w_brain_mask.nii.gz -o ${PTinHV_int[$p]}/sub-${HVs[$i]}_T1w_in_${PTs[$p]}_brain_mask.nii.gz -r ${PT_d_in[$p]}/sub-${PTs[$p]}_T1_nat_filled.nii.gz \
                    -t [${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_0GenericAffine.mat,1] -t ${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_1InverseWarp.nii.gz -n multilabel"

                    task_exec

                    task_in="fslmaths ${PTinHV_int[$p]}/sub-${HVs[$i]}_T1w_in_${PTs[$p]}_brain_mask.nii.gz -sub ${PT_intd[$p]}/sub-${PTs[$p]}_LM_dst.nii.gz ${PT_intd[$p]}/sub-${PTs[$p]}_LM_dst_inv.nii.gz"
                    
                    task_exec

                    sleep 5

                    task_in="mrcalc -force -quiet -nthreads ${ncpu} ${PT_intd[$p]}/sub-${PTs[$p]}_rT1w_brain.nii.gz \
                    ` mrstats -force -quiet -ignorezero -mask ${PT_intd[$p]}/sub-${PTs[$p]}_LM_dst_inv_BM.nii.gz -output mean ${PT_intd[$p]}/sub-${PTs[$p]}_rT1w_brain.nii.gz ` -div  \
                    ` mrstats -force -quiet -ignorezero -mask ${PTinHV_int[$p]}/sub-${HVs[$i]}_T1w_in_${PTs[$p]}_brain_mask.nii.gz -output mean ${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_InverseWarped.nii.gz ` \
                    -mult ${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_intmatched.nii.gz"

                    task_exec

                fi

                # exit 2

                if [[ ! -f "${PT_outad[$p]}/sub-${HVs[$i]}_in_${PTs[$p]}_SME_T1w.nii.gz" ]]; then

                    task_in="fslmaths ${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_intmatched.nii.gz \
                    -mul ${PT_intd[$p]}/sub-${PTs[$p]}_LM_dst.nii.gz ${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_lesion.nii.gz \
                    && fslmaths ${PTinHV_int[$p]}/sub-${HVs[$i]}_T1w_in_${PTs[$p]}_brain_mask.nii.gz -binv -mul ${PTinHV_int[$p]}/sub-${HVs[$i]}_T1w_in_${PTs[$p]}_Warped1.nii.gz \
                    ${PTinHV_int[$p]}/sub-${HVs[$i]}_T1w_in_${PTs[$p]}_skull.nii.gz \
                    && fslmaths ${PTinHV_int[$p]}/sub-${HVs[$i]}_T1w_in_${PTs[$p]}_Warped1.nii.gz -mul ${PT_intd[$p]}/sub-${PTs[$p]}_LM_dst_inv.nii.gz \
                    ${PTinHV_int[$p]}/sub-${HVs[$i]}_in_${PTs[$p]}_SP_T1w_brain_punched.nii.gz"

                    task_exec

                    task_in="fslmaths ${PTinHV_int[$p]}/sub-${HVs[$i]}_in_${PTs[$p]}_SP_T1w_brain_punched.nii.gz \
                    -add ${PTinHV_int[$p]}/sub-${PTs[$p]}_2${HVs[$i]}_lesion.nii.gz ${PTinHV_int[$p]}/sub-${HVs[$i]}_in_${PTs[$p]}_SP_T1w_brain.nii.gz \
                    && fslmaths ${PTinHV_int[$p]}/sub-${HVs[$i]}_in_${PTs[$p]}_SP_T1w_brain.nii.gz -add ${PTinHV_int[$p]}/sub-${HVs[$i]}_T1w_in_${PTs[$p]}_skull.nii.gz \
                    ${PT_outbd[$p]}/sub-${HVs[$i]}_in_${PTs[$p]}_SP_T1w.nii.gz \
                    && fslmaths ${PTinHV_int[$p]}/sub-${HVs[$i]}_T1w_in_${PTs[$p]}_Warped1.nii.gz -mas ${PTinHV_int[$p]}/sub-${HVs[$i]}_T1w_in_${PTs[$p]}_brain_mask.nii.gz \
                    -add ${PTinHV_int[$p]}/sub-${HVs[$i]}_T1w_in_${PTs[$p]}_skull.nii.gz ${PT_outad[$p]}/sub-${HVs[$i]}_in_${PTs[$p]}_SME_T1w.nii.gz"
                    
                    task_exec

                else

                    echo " sub-${HVs[$i]}_in_${PTs[$p]}_SME_T1w.nii.gz already generated, skipping to next one "

                fi

            fi

        done

    fi

done