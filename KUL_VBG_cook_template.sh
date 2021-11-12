#!/bin/bash

set -x

# Intro:
# - This script will cook a template for VBG specific to your population
# - This will ensure a better result from VBG
# - You only need to run this once for your study if you're using only 1 scanner
# - We use the whole head T1 WIs (any modality could work actually)
#   of two patients with opposing sided focal lesions 
#   (one with a left side lesion and one with right side lesion)
#   not causing mass effect nor midline shift
# - this makes it easier for VBG to generate a better fitting lesion fill
# 
# Instructions:
# 1- Edit input lines in this script to specify input images used for cooking your template
# 2- Place this script in the same directory as your inputs
# 3- Have KUL Neuroimaging tools in your path (e.g /KUL_NITs)
# 4- Have hd-bet installed and in your path (or simply edit as indicated below for using antsBrainExtraction instead of hd-bet)
# 5- After a successful run, place the cooked template images under KUL_NITs/atlasses/Templates/ 
# as VBG_T1_temp_brain.nii.gz and VBG_T1_temp.nii.gz

# Requirements:
# 1- T1 images from 2 patients with contralateral lesions without midline shift
# 2- To do: add option of using 2 brains with ipsilateral lesions

# Cook_template_4VBG.sh
# v=0.1

# AR, SS 14/04/2020

# Misc. vars
# change me to suit your CPU :)

ncpu=24

d=$(date "+%Y-%m-%d_%H-%M-%S");

cwd=$(pwd)

wkdir="${cwd}/VBG_cook_temp"

outdir="${cwd}/VBG_cooked_temp"

###### CHANGE ME #### CHANGE ME
# T1_R is the one with a healthy right side
# T1_L is the one with a healthy left side

# T1_R="${cwd}/sub-001_T1w.nii.gz"

# T1_L="${cwd}/sub-001_T1w_F11.nii.gz"

T1_R="/NI_apps/KUL_VBG/share/Test_data/Pats/sub-PT010/sub-PT010_T1w_defaced.nii.gz"

T1_L="/NI_apps/KUL_VBG/share/Test_data/Pats/sub-PT012/sub-PT012_T1w_defaced.nii.gz"

###### CHANGE ME #### CHANGE ME

# make dirs

mkdir -p ${wkdir}

mkdir -p ${outdir}

# prep log

prep_log="${wkdir}/prep_log_${d}.txt";

if [[ ! -f ${prep_log} ]] ; then

    touch ${prep_log}

else

    echo "${prep_log} already created"

fi

# Temp. vars

temps_dir=($(which KUL_VBG.sh | rev | cut -d"/" -f2- | rev))

MNI_T1="${temps_dir}/atlasses/Templates/HR_T1_MNI.nii.gz"

MNI_T1_brain="${temps_dir}/atlasses/Templates/HR_T1_MNI_brain.nii.gz"

MNI_brain_mask="${temps_dir}/atlasses/Templates/HR_T1_MNI_brain_mask.nii.gz"

MNI_r="${temps_dir}/atlasses/Templates/Rt_hemi_mask.nii.gz"

MNI_l="${temps_dir}/atlasses/Templates/Lt_hemi_mask.nii.gz"

# Workflow vars

# T1_Rf="${wkdir}/sub-PT024_T1w_flipped.nii.gz"

HDBET_L="${wkdir}/T1_L_HDBET"

HDBET_R="${wkdir}/T1_R_HDBET"

T1_L_brain_2MNI1="${wkdir}/T1_L_HDBET_2MNI1"

T1_R_brain_2MNI1="${wkdir}/T1_R_HDBET_2MNI1"

T1_L_inMNI="${wkdir}/T1_L_inMNI1.nii.gz"

T1_L_brain_inMNI="${wkdir}/T1_L_HDBET_2MNI1_Warped.nii.gz"

T1_R_brain_inMNI="${wkdir}/T1_R_HDBET_2MNI1_Warped.nii.gz"

MNI2_in_T1="${wkdir}/T1_brain_inMNI2_InverseWarped.nii.gz"

MNI_lw="${wkdir}/MNI_L_insubjT1_inMNI1.nii.gz"

MNI_lwr="${wkdir}/MNI_L_insubjT1_inMNI1r.nii.gz"

MNI_rwr="${wkdir}/MNI_R_insubjT1_inMNI1r.nii.gz"

T1_L_skull="${wkdir}/T1_L_inMNI1_skull.nii.gz"

T1_L_skull_mask="${wkdir}/T1_L_inMNI1_skull_mask.nii.gz"

T1_L_brain_2MNI2="${wkdir}/T1_L_brain_2MNI2"

T1_R_brain_2MNI2="${wkdir}/T1_R_brain_2MNI2"

MNI_skull_mask="${wkdir}/MNI_skull_mask.nii.gz"

MNI_skull_HM="${wkdir}/MNI_skull_HM2_T1_L.nii.gz"

T1_L_brain_inMNI2_HM2R_b="${wkdir}/T1_L_brain_HM2_T1_R_blurry.nii.gz"

T1_L_brain_inMNI2_HM2R="${wkdir}/T1_L_brain_HM2_T1_R.nii.gz"

ARZ_brain="${wkdir}/ARZ_T1_brain.nii.gz"

ARZ_T1="${wkdir}/ARZ_T1.nii.gz"

cooked_template_brain="${outdir}/VBG_T1_temp_brain.nii.gz"

cooked_template_T1="${outdir}/VBG_T1_temp.nii.gz"

# task exec function

function task_exec {

    echo "  " >> ${prep_log} 
    
    echo ${task_in} >> ${prep_log} 

    echo " Started @ $(date "+%Y-%m-%d_%H-%M-%S")" >> ${prep_log} 

    eval ${task_in} >> ${prep_log} 2>&1 &

    echo " pid = $! basicPID = $$ " >> ${prep_log}

    wait ${pid}

    sleep 5

    if [ $? -eq 0 ]; then
        echo Success >> ${prep_log}
    else
        echo Fail >> ${prep_log}

        exit 1
    fi

    echo " Finished @  $(date "+%Y-%m-%d_%H-%M-%S")" >> ${prep_log} 

    echo "  " >> ${prep_log} 

    unset task_in

}

# please cite us

echo " You are using the cook template script, part of the VBG package, please cite the following paper in your work: \
Virtual brain grafting: Enabling whole brain parcellation in the presence of large lesions \
Ahmed M. Radwan, Louise Emsell, Jeroen Blommaert, Andrey Zhylka, Silvia Kovacs, Tom Theys, Nico Sollmann, Patrick Dupont, Stefan Sunaert \
medRxiv 2020.09.30.20204701; doi: https://doi.org/10.1016/j.neuroimage.2021.117731" | tee -a ${prep_log}

# 1- HDBET both T1s and use warps to get to MNI

# task_in="antsBrainExtraction.sh -d 3 -a ${T1_L} -e ${MNI_T1} -m ${MNI_brain_mask} -u 1 -o ${HDBET_L}_"
# if no hd-bet then change the ${HDBET_L}, ${HD_BET_R}, and consequent command inputs to suit your brain extraction tool.


# task_exec

task_in="hd-bet -i ${T1_L} -o ${HDBET_L}_brain"

task_exec

task_in="hd-bet -i ${T1_R} -o ${HDBET_R}_brain"

task_exec

#

task_in="antsRegistrationSyN.sh -n ${ncpu} -d 3 -f ${MNI_T1_brain} -m ${HDBET_L}_brain.nii.gz -x ${MNI_brain_mask},${HDBET_L}_brain_mask.nii.gz -t s -n ${ncpu} -o ${T1_L_brain_2MNI1}_"

task_exec &

task_in="antsRegistrationSyN.sh -n ${ncpu} -d 3 -f ${MNI_T1_brain} -m ${HDBET_R}_brain.nii.gz -x ${MNI_brain_mask},${HDBET_R}_brain_mask.nii.gz -t s -n ${ncpu} -o ${T1_R_brain_2MNI1}_"

task_exec

task_in="WarpImageMultiTransform 3 ${T1_L} ${T1_L_inMNI} -R ${MNI_T1}  ${T1_L_brain_2MNI1}_1Warp.nii.gz ${T1_L_brain_2MNI1}_0GenericAffine.mat"

task_exec


task_in="mrcalc -quiet -force -nthreads ${ncpu} ${T1_L_brain_inMNI} -neg 0 -ge ${T1_L_inMNI} -mult ${T1_L_skull} && mrcalc -quiet -force -nthreads ${ncpu} ${T1_L_skull} 0 -gt ${T1_L_skull_mask}"

task_exec

# 2- antsRegSyN 1 stage to MNI again

task_in="antsRegistrationSyN.sh -d 3 -f ${MNI_T1_brain} -m ${T1_L_brain_inMNI} -x ${MNI_brain_mask},${MNI_T1_brain_mask} -t so -n ${ncpu} -o ${T1_L_brain_2MNI2}_"

task_exec &

task_in="antsRegistrationSyN.sh -d 3 -f ${MNI_T1_brain} -m ${T1_R_brain_inMNI} -x ${MNI_brain_mask},${MNI_T1_brain_mask} -t so -n ${ncpu} -o ${T1_R_brain_2MNI2}_"

task_exec

task_in="ImageMath 3 ${T1_R_brain_2MNI2}_sharp_Warped.nii.gz Sharpen ${T1_R_brain_2MNI2}_Warped.nii.gz"

task_exec

# 3- apply inverse syn warp to MNI R/L masks
# this step should probably use the prior from only one brain

task_in="WarpImageMultiTransform 3 ${MNI_l} ${MNI_lw} -R ${T1_L_brain_inMNI} -i ${T1_L_brain_2MNI2}_0GenericAffine.mat ${T1_L_brain_2MNI2}_1InverseWarp.nii.gz"

task_exec

task_in="mrcalc -quiet -force -nthreads ${ncpu} ${MNI_lw} ${MNI_brain_mask} -mult 0.1 -gt ${MNI_lwr} && mrcalc -quiet -force -nthreads ${ncpu} ${MNI_lwr} -neg 0 -ge ${MNI_brain_mask} -mult ${MNI_rwr}"

task_exec

# task_in="WarpImageMultiTransform 3 ${MNI_r} ${MNI_rw} -R ${T1_R_brain_inMNI} -i ${T1_R_brain_2MNI2}_1InverseWarp.nii.gz ${T1_R_brain_2MNI2}_0GenericAffine.mat"

# task_exec

# 4- Make MNI skull and hist match to subject skull

task_in="mrcalc -quiet -force -nthreads ${ncpu} ${MNI_brain_mask} -neg 0 -ge ${MNI_T1} -mult 0 -gt ${MNI_skull_mask}" 

task_exec

task_in="mrcalc -quiet -force -nthreads ${ncpu} ${MNI_brain_mask} -neg 0 -ge ${MNI_T1} -mult - | mrhistmatch linear - ${T1_L_skull} ${MNI_skull_HM} -nthreads ${ncpu} \
-mask_target ${T1_L_skull_mask} -mask_input ${MNI_skull_mask} -quiet"

task_exec

# 4- Hist match brains in MNI

task_in="mrhistmatch -quiet -bins 2048 -force -nthreads ${ncpu} -mask_input ${MNI_brain_mask} -mask_target ${MNI_brain_mask} nonlinear ${T1_L_brain_2MNI2}_Warped.nii.gz \
${T1_R_brain_2MNI2}_Warped.nii.gz ${T1_L_brain_inMNI2_HM2R_b} && ImageMath 3 ${T1_L_brain_inMNI2_HM2R} Sharpen ${T1_L_brain_inMNI2_HM2R_b}"

task_exec

# 5- Mask and stitch healthy hemis to cook template ims

task_in="mrcalc -force -nthreads ${ncpu} ${T1_L_brain_inMNI2_HM2R} ${MNI_lwr} -mult `mrcalc -quiet -force -nthreads ${ncpu} ${T1_R_brain_2MNI2}_sharp_Warped.nii.gz ${MNI_rwr} -mult - ` \
-add ${ARZ_brain} && mrcalc -force -nthreads ${ncpu} ${ARZ_brain} ${MNI_brain_mask} -mult ${MNI_skull_HM} -add ${ARZ_T1}"

task_exec

task_in="WarpImageMultiTransform 3 ${ARZ_T1} ${cooked_template_T1} -R ${MNI_T1_brain} -i ${T1_L_brain_2MNI2}_0GenericAffine.mat ${T1_L_brain_2MNI2}_1InverseWarp.nii.gz && WarpImageMultiTransform 3 \
${ARZ_brain} ${cooked_template_brain} -R ${MNI_T1_brain} -i ${T1_L_brain_2MNI2}_0GenericAffine.mat ${T1_L_brain_2MNI2}_1InverseWarp.nii.gz"

task_exec

# done

