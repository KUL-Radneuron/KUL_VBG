#!/bin/bash

# set -x

# This script will grab HV T1s, reorient them and extract brain with Bet
# it will then reg a patient brain to HV brain
# and apply warp to lesion mask with --NN interpolation
# this is preparing for running FS for all these HVs

# inputs needed:
# from patients:
# 1- filled T1 brain 2- T1 brain with lesion 3- T1 brain mask 4- Lmask_reori_bin
# from healthy volunteers:
# 1- whole head T1

# to do:
# 1- define vars and start for loop for HVs
# 2- Run reorient2std and HD-BET for HVs
# 3- nested for loop for PTs
# 4- Reg each PT to each HV
# 5- apply inv warps to bring HV head to PT 
# 6- generate Lmask binv
# 7- intensity match PT brain with lesion to HV brain
# 8- transplant the lesion into HV head in PT space

# dirs needed:
# 1- work dir
# 2- input for PT images, subfolder per patient is a good idea
# should contain T1 brain with lesion, filled, lesion mask reori bin, PT brain mask
# 3- input for HVs, can be a mess
# 4- input for masks
# 5- processing dir should contain a- registration dir, with per sub reg dirs
# b- lesion_gp dirs ?

ncpu="36"

wd=$(pwd)

HV_dir="${wd}/HVs_NII"

Pats_d="${wd}/Pat_T1s"
# this should contain the filled brains, and brains with lesions of the patients

# so we have one parent folder, one folder for GT, one for test data and one for processing
val_dir="${wd}/KUL_VBG_validation"

wrk_dir="${val_dir}/Proc_dir"

regs_dir="${wrk_dir}/Pats_2HVs_reg_dir"

work_dir_masks="${regs_dir}/Lesion_masks_proc"

test_dir="${val_dir}/VBG_test"

test_dir_T1s="${val_dir}/VBG_test/T1s"

# test_dir_masks="${val_dir}/VBG_test/Lesion_masks"

GT_dir="${val_dir}/VBG_GT"

# dice_dir

# make your dirs

mkdir -p ${wrk_dir}

mkdir -p ${val_dir}

mkdir -p ${test_dir}

mkdir -p ${test_dir_T1s}

mkdir -p ${regs_dir}

mkdir -p ${work_dir_masks}

mkdir -p ${test_dir_masks}

mkdir -p ${GT_dir}

HVs=("HV001"  "HV002"  "HV003"  "HV004"  "HV005"  "HV006"  "HV007"  \
"HV008" "HV009"  "HV010");

PTs=("PT051" "PT010" "PT014" "PT043" "PT034" "PT020" \
"PT024" "PT052" "PT012" "PT021");

# PTs=("PT020" "PT024" "PT052" "PT012" "PT021")

HV_T1s=("20181025071243_CS_T1_TFE_0.9mm_501.nii.gz"  "20181119180029_CS_T1_TFE_0.9mm_301.nii.gz"  \
"20181122100602_WIP_CS_T1_TFE_0.9mm_301.nii.gz"  "20181122111817_WIP_CS_T1_TFE_0.9mm_201.nii.gz"  \
"20181122181018_CS_T1_TFE_0.9mm_301.nii.gz"  "20181126091638_WIP_CS_T1_TFE_0.9mm_201.nii.gz"  \
"20181203182250_CS_T1_TFE_0.9mm_201.nii.gz"  "20181204123520_WIP_CS_T1_TFE_0.9mm_201.nii.gz"  \
"20181219172827_CS_T1_TFE_0.9mm_201.nii.gz"  "20190108191325_CS_T1_TFE_0.9mm_401.nii.gz");

d=$(date "+%Y-%m-%d_%H-%M-%S");

sub_log="${wrk_dir}/trial_${d}.txt"

touch "${sub_log}"

function task_exec {

    echo "-------------------------------------------------------------" >> ${sub_log}

    echo ${task_in} >> ${sub_log}

    echo " Started @ $(date "+%Y-%m-%d_%H-%M-%S")" >> ${sub_log}

    eval ${task_in} >> ${sub_log} 2>&1 &

    echo " pid = $! basicPID = $$ " >> ${sub_log}

    wait ${pid}

    sleep 5

    if [ $? -eq 0 ]; then

        echo Success >> ${sub_log}

    else

        echo Fail >> ${sub_log}

        exit 1

    fi

    echo " Finished @ $(date "+%Y-%m-%d_%H-%M-%S")" >> ${sub_log}

    echo "-------------------------------------------------------------" >> ${sub_log}
    
    echo "" >> ${sub_log}

    unset task_in


}

####

# make a bin and binv Lmask # should apply a slight smoothing also

# N=10

for pt in ${!PTs[@]}; do

    PT_masks="${work_dir_masks}/sub-${PTs[$pt]}_masks"

    mkdir -p ${PT_masks}

    # ((pa=pa%N)); ((pa++==0)) && wait

    if [[ ! -f "${PT_masks}/sub-${PTs[$pt]}_WLE_AR_nat_dilms2_inv_BM.nii.gz" ]]; then 

        task_in="fslmaths ${Pats_d}/sub-${PTs[$pt]}/sub-${PTs[$pt]}_WLE_AR_nat.nii.gz -dilM -s 2 -thr 0.2 -save \
        ${PT_masks}/sub-${PTs[$pt]}_WLE_AR_nat_dilms2.nii.gz -binv -mul ${Pats_d}/sub-${PTs[$pt]}/sub-${PTs[$pt]}_Gs_ARn_Brain_clean_mask.nii.gz \
        ${PT_masks}/sub-${PTs[$pt]}_WLE_AR_nat_dilms2_inv_BM.nii.gz"

        task_exec &

    else

        echo " ${PT_masks}/sub-${PTs[$pt]}_WLE_AR_nat_dilms2_inv.nii.gz already done"

    fi

done

wait


### deal with the HVs

# k=3

for i in ${!HVs[@]}; do

    # ((pb=pb%k)); ((pa++==0)) && wait

    HV_wrk="${regs_dir}/sub-${HVs[$i]}_work"

    mkdir -p ${HV_wrk}

    echo "Now working on HVs BETs"

    if [[ ! -f "${HV_wrk}/sub-${HVs[$i]}_T1_brain.nii.gz" ]]; then

        task_in="fslreorient2std ${HV_dir}/${HV_T1s[$i]} ${HV_wrk}/sub-${HVs[$i]}_T1_reori.nii.gz"

        task_exec
    
        task_in="hd-bet -i ${HV_wrk}/sub-${HVs[$i]}_T1_reori.nii.gz -o ${HV_wrk}/sub-${HVs[$i]}_T1_brain"

        task_exec

    else

        echo " ${HV_wrk}/sub-${HVs[$i]}_T1_brain already done"

    fi

    # here we use the filled patient images to achieve best spatial match to the healthy volunteers
    # we use default parameters for ANTsregSyN
    # the warps will be applied to bring the lesion masks from patient space to HV space
    # we will continue working only with the subject images in their original space
    # after everything is finished we can bring the aparc+aseg_minL.nii.gz to patient space 
    # we treat each patient's space as a reference space for the lesion type

    for p in ${!PTs[@]}; do

        # ((pa=pa%N)); ((pa++==0)) && wait

        # better make different folders for each PT gp. of HVs
        # forget about looking at deformed and non-deformed, we only have 1 ground truth - aparc+aseg from deformed nonlesioned HVs

        PT_d_in="${Pats_d}/sub-${PTs[$p]}"

        PT_d_w="${regs_dir}/sub-${PTs[$p]}_HVs"

        PT_masks="${work_dir_masks}/sub-${PTs[$p]}_masks"

        mkdir -p ${PT_d_w}

        mkdir -p ${PT_masks}

        # ants on steroids

        if [[ ! -f "${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_Warped.nii.gz" ]]; then

            # task_in="export ITK_GLOBAL_DEFAULT_NUMBER_OF_THRclearEADS=${ncpu} ; antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 -u 1 \
            # --output [ ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_,${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_Warped.nii.gz,${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_InverseWarped.nii.gz ] \
            # --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [ 0.005,0.995 ] \
            # -x [ ${HV_wrk}/sub-${HVs[$i]}_T1_brain_mask.nii.gz,${PT_d_in}/sub-${PTs[$p]}_Gs_ARn_Brain_clean_mask.nii.gz, NULL ] \
            # --initial-moving-transform [ ${HV_wrk}/sub-${HVs[$i]}_T1_brain.nii.gz,${PT_d_in}/sub-${PTs[$p]}_T1_nat_filled.nii.gz,1 ] --transform Rigid[ 0.1 ] \
            # --metric MI[ ${HV_wrk}/sub-${HVs[$i]}_T1_brain.nii.gz,${PT_d_in}/sub-${PTs[$p]}_T1_nat_filled.nii.gz,1,32,Regular,0.25 ] \
            # --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[ 0.1 ] \
            # --metric MI[ ${HV_wrk}/sub-${HVs[$i]}_T1_brain.nii.gz,${PT_d_in}/sub-${PTs[$p]}_T1_nat_filled.nii.gz,1,64,Regular,0.5 ] \
            # --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform SyN[ 0.1,3,0 ] \
            # --metric CC[ ${HV_wrk}/sub-${HVs[$i]}_T1_brain.nii.gz,${PT_d_in}/sub-${PTs[$p]}_T1_nat_filled.nii.gz,1,4 ] \
            # --convergence [ 200x100x75x25,1e-8,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --verbose 1"

            echo "export ITK_GLOBAL_DEFAULT_NUMBER_OF_THRclearEADS=${ncpu} ; antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 -u 1 \
            --output [ ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_,${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_Warped.nii.gz,${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_InverseWarped.nii.gz ] \
            --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [ 0.005,0.995 ] \
            -x [ ${HV_wrk}/sub-${HVs[$i]}_T1_brain_mask.nii.gz,${PT_d_in}/sub-${PTs[$p]}_Gs_ARn_Brain_clean_mask.nii.gz, NULL ] \
            --initial-moving-transform [ ${HV_wrk}/sub-${HVs[$i]}_T1_brain.nii.gz,${PT_d_in}/sub-${PTs[$p]}_T1_nat_filled.nii.gz,1 ] --transform Rigid[ 0.1 ] \
            --metric MI[ ${HV_wrk}/sub-${HVs[$i]}_T1_brain.nii.gz,${PT_d_in}/sub-${PTs[$p]}_T1_nat_filled.nii.gz,1,32,Regular,0.25 ] \
            --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[ 0.1 ] \
            --metric MI[ ${HV_wrk}/sub-${HVs[$i]}_T1_brain.nii.gz,${PT_d_in}/sub-${PTs[$p]}_T1_nat_filled.nii.gz,1,64,Regular,0.5 ] \
            --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform SyN[ 0.1,3,0 ] \
            --metric CC[ ${HV_wrk}/sub-${HVs[$i]}_T1_brain.nii.gz,${PT_d_in}/sub-${PTs[$p]}_T1_nat_filled.nii.gz,1,4 ] \
            --convergence [ 200x100x75x25,1e-8,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --verbose 1"  >> ${sub_log}

            # task_exec

        else

            echo " sub-${PTs[$p]}_2${HVs[$i]}_Warped.nii.gz already generated, skipping to next one"

        fi

            #  --verbose 1

            # 1 task in to rule them all

            # task_in="WarpImageMultiTransform 3 ${HV_wrk}/sub-${HVs[$i]}_T1_reori.nii.gz ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}.nii.gz  \
            # -R ${PT_d_in}/sub-${PTs[$p]}_T1_nat_filled.nii.gz -i ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_0GenericAffine.mat \
            # ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_1InverseWarp.nii.gz && \
            # \
            # WarpImageMultiTransform 3 ${HV_wrk}/sub-${HVs[$i]}_T1_brain_mask.nii.gz \
            # ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask.nii.gz -R ${PT_d_in}/sub-${PTs[$p]}_T1_nat_filled.nii.gz \
            # -i ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_0GenericAffine.mat ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_1InverseWarp.nii.gz && \
            # \
            # fslmaths ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask.nii.gz -thr 0.5 -bin -save ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask_bin.nii.gz \
            # -sub ${PT_masks}/sub-${PTs[$p]}_WLE_AR_nat_dilms2.nii.gz ${PT_masks}/sub-${PTs[$p]}_WLE_AR_nat_dilms2_inv.nii.gz && \
            # \
            # mrcalc -force -nthreads ${ncpu} ${PT_d_in}/sub-${PTs[$p]}_Gs_ARn_Brain_clean.nii.gz \
            # ` mrstats -ignorezero -quiet -mask ${PT_masks}/sub-${PTs[$p]}_WLE_AR_nat_dilms2_inv_BM.nii.gz -output mean ${PT_d_in}/sub-${PTs[$p]}_Gs_ARn_Brain_clean.nii.gz ` -div \
            # ` mrstats -ignorezero -quiet -mask ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask_bin.nii.gz -output mean ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_InverseWarped.nii.gz ` -mult \
            # ${PT_d_w}/sub-${PTs[$p]}_T1_in_${HVs[$i]}_intmatch.nii.gz && \
            # \
            # fslmaths ${PT_d_w}/sub-${PTs[$p]}_T1_in_${HVs[$i]}_intmatch.nii.gz \
            # -mul ${PT_masks}/sub-${PTs[$p]}_WLE_AR_nat_dilms2.nii.gz ${PT_masks}/sub-${PTs[$p]}_T1_in_${HVs[$i]}_lesion.nii.gz && \
            # \
            # fslmaths ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask.nii.gz -binv -mul ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}.nii.gz -save \
            # ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_skull.nii.gz \
            # \
            # -restart ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}.nii.gz -mul \
            # ${PT_masks}/sub-${PTs[$p]}_WLE_AR_nat_dilms2_inv.nii.gz -add ${PT_masks}/sub-${PTs[$p]}_T1_in_${HVs[$i]}_lesion.nii.gz \
            # -mas ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask.nii.gz -add ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_skull.nii.gz -save \
            # ${test_dir_T1s}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_lesioned.nii.gz \
            # \
            # -restart ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}.nii.gz -mul \
            # ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask.nii.gz -add ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_skull.nii.gz \
            # ${GT_dir}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_GT.nii.gz"

            # task_exec &

            # task_exec

            echo "WarpImageMultiTransform 3 ${HV_wrk}/sub-${HVs[$i]}_T1_reori.nii.gz ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}.nii.gz  \
            -R ${PT_d_in}/sub-${PTs[$p]}_T1_nat_filled.nii.gz -i ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_0GenericAffine.mat \
            ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_1InverseWarp.nii.gz"  >> ${sub_log}
            
            echo "WarpImageMultiTransform 3 ${HV_wrk}/sub-${HVs[$i]}_T1_brain_mask.nii.gz \
            ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask.nii.gz -R ${PT_d_in}/sub-${PTs[$p]}_T1_nat_filled.nii.gz \
            -i ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_0GenericAffine.mat ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_1InverseWarp.nii.gz"  >> ${sub_log}
            
            echo "fslmaths ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask.nii.gz -thr 0.5 -bin -save ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask_bin.nii.gz \
            -sub ${PT_masks}/sub-${PTs[$p]}_WLE_AR_nat_dilms2.nii.gz ${PT_masks}/sub-${PTs[$p]}_WLE_AR_nat_dilms2_inv.nii.gz"  >> ${sub_log}
            
            echo "mrcalc -force -nthreads ${ncpu} ${PT_d_in}/sub-${PTs[$p]}_Gs_ARn_Brain_clean.nii.gz \
            ` mrstats -ignorezero -quiet -mask ${PT_masks}/sub-${PTs[$p]}_WLE_AR_nat_dilms2_inv_BM.nii.gz -output mean ${PT_d_in}/sub-${PTs[$p]}_Gs_ARn_Brain_clean.nii.gz ` -div \
            ` mrstats -ignorezero -quiet -mask ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask_bin.nii.gz -output mean ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_InverseWarped.nii.gz ` -mult \
            ${PT_d_w}/sub-${PTs[$p]}_T1_in_${HVs[$i]}_intmatch.nii.gz"  >> ${sub_log}

            echo " mrstats -ignorezero -quiet -mask ${PT_masks}/sub-${PTs[$p]}_WLE_AR_nat_dilms2_inv_BM.nii.gz -output mean ${PT_d_in}/sub-${PTs[$p]}_Gs_ARn_Brain_clean.nii.gz "  >> ${sub_log}
 
            echo "mrstats -ignorezero -quiet -mask ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask_bin.nii.gz -output mean ${PT_d_w}/sub-${PTs[$p]}_2${HVs[$i]}_InverseWarped.nii.gz"  >> ${sub_log}
            
            echo "fslmaths ${PT_d_w}/sub-${PTs[$p]}_T1_in_${HVs[$i]}_intmatch.nii.gz \
            -mul ${PT_masks}/sub-${PTs[$p]}_WLE_AR_nat_dilms2.nii.gz ${PT_masks}/sub-${PTs[$p]}_T1_in_${HVs[$i]}_lesion.nii.gz"  >> ${sub_log}
            
            echo "fslmaths ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask.nii.gz -binv -mul ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}.nii.gz -save \
            ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_skull.nii.gz \
            -restart ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}.nii.gz -mul \
            ${PT_masks}/sub-${PTs[$p]}_WLE_AR_nat_dilms2_inv.nii.gz -add ${PT_masks}/sub-${PTs[$p]}_T1_in_${HVs[$i]}_lesion.nii.gz \
            -mas ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask.nii.gz -add ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_skull.nii.gz -save \
            ${test_dir_T1s}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_lesioned.nii.gz \
            \
            -restart ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}.nii.gz -mul \
            ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_brain_mask.nii.gz -add ${PT_d_w}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_skull.nii.gz \
            ${GT_dir}/sub-${HVs[$i]}_T1_in_${PTs[$p]}_GT.nii.gz"  >> ${sub_log}


        else

            echo " sub-${PTs[$p]}_2${HVs[$i]}_Warped.nii.gz already generated, skipping to next one"

        fi

        # # now we are ready for VBG 
        # # first we run VBG_BETonly
        # # then VBG for all in parallel (this we do locally)
        # # HVs_inPT T1 (deformed nonlesioned) and regular FS on the lesioned HVs without VBG (to record success v failure)
        # # do the last one in a separate script

    done

    # these commands might come in handy later on

    # task_in="antsRegistrationSyN.sh -d 3 -f /NI_apps/FSL/data/standard/MNI152_T1_1mm_brain.nii.gz -m ${work_dir}/${HVs[$i]}_T1_brain.nii.gz \
    # -x /NI_apps/FSL/data/standard/MNI152_T1_1mm_brain_mask.nii.gz,${work_dir}/${HVs[$i]}_T1_brain_mask.nii.gz -t a -n ${ncpu} \
    # -o ${work_dir}/${HVs[$i]}_T1_brain_2MNIa_"

    # task_exec

    # task_in="WarpImageMultiTransform 3 ${work_dir}/${HVs[$i]}_T1_brain_mask.nii.gz ${work_dir}/${HVs[$i]}_T1_bm_inMNIa.nii.gz  \
    # -R /NI_apps/FSL/data/standard/MNI152_T1_1mm_brain.nii.gz ${work_dir}/${HVs[$i]}_T1_brain_2MNIa_0GenericAffine.mat"

    # task_exec

    # labelconvert -force -nthreads 2 ${wd}/lesion_wf/output_LWF/sub-${Gs_sM[$i]}_Gs_ARt_aparc+aseg_minL.nii.gz \
    # /NI_apps/freesurfer/FreeSurferColorLUT.txt /NI_apps/mrtrix3/share/mrtrix3/labelconvert/fs_default.txt \
    # ${wd}/lesion_wf/output_LWF/sub-${Gs_sM[$i]}_Gs_ARt_aparc+aseg_minL_LC.nii.gz 2>&1 > ${run_log_t} &


done
