#!/bin/bash

# set -x

# Ahmed Radwan ahmed.radwan@kuleuven.be
# Stefan Sunaert  stefan.sunaert@kuleuven.be

# this script is in dev for S61759

#####################################


v="0.51_13052021_beta"

# This script is meant to allow a decent recon-all/antsMALF output in the presence of a large brain lesion 
# The main idea is to replace the lesion with a hole and fill the hole with information from the a synthetic image
# this maintains subject specificity and diseased hemisphere information but replaces lesioned tissue with sham brain 
# to do:
# (1) Add option for FreeSurfer or FastSurfer ;)

# ----------------------------------- MAIN --------------------------------------------- 
# this script uses "preprocessing control", i.e. if some steps are already processed it will skip these

kul_lesion_dir=`dirname "$0"`
script=`basename "$0"`
# source $kul_main_dir/KUL_main_functions.sh
cwd=($(pwd))

# FUNCTIONS --------------

# function Usage
function Usage {

cat <<USAGE

`basename $0` preps structural images with lesions and runs recon-all.

Usage examples:

    `basename $0` -S Subject_ID -l /fullpath/lesion_mask.nii.gz -z T1 -b -B (1 or 2)
  
    or
  
    `basename $0` -S Subject_ID -a /fullpath/Subject_ID_T1w.nii.gz -l /fullpath/lesion_mask.nii.gz -z T1 -B (1 or 2) -n 10 -v 1

Purpose:

    The purpose of this workflow is to generate a lesion filled image, with healthy looking synthetic tissue in place of the lesion
    Essentially excising the lesion and grafting over the brain tissue defect in the MR image space

How to use:

    - You need to use the cook_template_4VBG script once for your study - if you have only 1 scanner
    - cook_template_4VBG requires two brains with unilateral lesions on opposing sides
    - it is meant to facilitate the grafting process and minimize intensity differences
    - You need a high resolution T1 WI and a lesion mask in the same space for VBG to run
    - If you end up with an empty image, it is possible you have a mismatch between the T1 and lesion mask


Required arguments:

    -S:  BIDS subject/participant name (anonymised name of the subject without the "sub-" prefix)
    -b:  if data is in BIDS
    -l:  full path and file name to lesion mask file per session
    -z:  space of the lesion mask used (only T1 supported in this version)
    -a:  Input precontrast T1WIs


Optional arguments:

    -s:  session (of the participant)
    -t:  Use the VBG template to derive the fill patch (if used, template tissue is used alongside native tissue to create the donor brain)
    -E:  Treat as an extra-axial lesion (skip VBG bulk, fill lesion patch with 0s, run FS and subsequent steps)
    -B:  specify brain extraction method (1 = HD-BET, 2 = ANTs-BET), if not set ANTs-BET will be used by default
    -P:  Run parcellation (1 = FreeSurfer recon-all, 2 = FastSurfer)
    -p:  In case of pediatric patients - use pediatric template (NKI_under_10 in MNI)
    -m:  full path to intermediate output dir
    -o:  full path to output dir (if not set reverts to default output ./VBG_output)
    -n:  number of cpu for parallelisation (default is 6)
    -v:  show output from mrtrix commands
    -h:  prints help menu

Notes: 

    - You can use -b and the script will find your BIDS files automatically
    - If your data is not in BIDS, then use -a without -b
    - This version is for validation only.
    - In case of trouble with HD-BET see lines (1124 - 1200)



USAGE

    exit 1
}


# CHECK COMMAND LINE OPTIONS -------------
# 
# Set defaults
# this works for ANTsX scripts and FS
ncpu=8


# Set required options
S_flag=0
b_flag=0
s_flag=0
l_flag=0
l_spaceflag=0
t1_flag=0
t_flag=0
o_flag=0
m_flag=0
n_flag=0
P_flag=0
E_flag=0
p_flag=0

if [ "$#" -lt 1 ]; then
    Usage >&2
    exit 1

else

    while getopts "S:a:l:z:s:o:m:n:B:P:bvhtEp" OPT; do

        case $OPT in
        S) #subject
            S_flag=1
            subj=$OPTARG
        ;;
        b) #BIDS or not ?
            bids_flag=1
        ;;
        a) #T1 WIs
            t1_flag=1
			t1_orig=$OPTARG
        ;;
        s) #session
            s_flag=1
            ses=$OPTARG
        ;;
        B) #session
            BET_flag=1
            BET_m=$OPTARG
        ;;
        l) #lesion_mask
            l_flag=1
            L_mask=$OPTARG
		;;
	    z) #lesion_mask
	        l_spaceflag=1
	        L_mask_space=$OPTARG	
	    ;;
	    m) #intermediate output
			m_flag=1
			wf_dir=$OPTARG		
        ;;
	    o) #output
			o_flag=1
			out_dir=$OPTARG		
        ;;
        t) #template flag
			t_flag=1	
        ;;
        P) #Parcellation flag
			P_flag=1
            parc_F=$OPTARG
        ;;
        E) #Extra-axial flag
			E_flag=1	
        ;;
        p) #Extra-axial flag
			p_flag=1	
        ;;
        n) #parallel
			n_flag=1
            ncpu=$OPTARG
        ;;
        v) #verbose
            silent=0
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

lesion_wf="${cwd}/VBG_out"

# output

if [[ "$o_flag" -eq 1 ]]; then
	
    output_m="${out_dir}"

    output_d="${output_m}/output_VBG/sub-${subj}${ses_long}"

else

	output_d="${lesion_wf}/output_VBG/sub-${subj}${ses_long}"

fi

# intermediate folder

if [[ "$m_flag" -eq 1 ]]; then

	preproc_m="${wf_dir}"

    preproc="${preproc_m}/proc_VBG/sub-${subj}${ses_long}"

else

	preproc="${lesion_wf}/proc_VBG/sub-${subj}${ses_long}"

fi


# timestamp
start_t=$(date +%s)

FSLPARALLEL=$ncpu; export FSLPARALLEL
OMP_NUM_THREADS=$ncpu; export OMP_NUM_THREADS

d=$(date "+%Y-%m-%d_%H-%M-%S");

# handle the dirs

cd $cwd

long_bids_subj="${search_sessions}"

echo $long_bids_subj

bids_subj=${long_bids_subj%anat}

echo $bids_subj

######

# check for required inputs and define your workflow accordingly

srch_Lmask_str=($(basename ${L_mask}))
srch_Lmask_dir=($(dirname ${L_mask}))
srch_Lmask_o=($(find ${srch_Lmask_dir} -type f | grep  ${srch_Lmask_str}))

if [[ $S_flag -eq 0 ]] || [[ $l_flag -eq 0 ]] || [[ $l_spaceflag -eq 0 ]]; then
	
    echo
    echo "Inputs -S -l -z must be set." >&2
    echo
    exit 2
	
else

    if [[ -z "${srch_Lmask_o}" ]]; then

        echo
        echo " Incorrect Lesion mask, please check the file path and name "
        echo
        exit 2

    else
	
	    echo "Inputs are -S  ${subj}  -l  ${L_mask}  -z  ${L_mask_space}"


    fi
	
fi

	
if [[ "$bids_flag" -eq 1 ]] && [[ "$s_flag" -eq 0 ]]; then
		
	# bids flag defined but not session flag
    search_sessions=($(find ${cwd}/BIDS/sub-${subj} -type d | grep anat));
	num_sessions=${#search_sessions[@]};
	ses_long="";
	
	if [[ "$num_sessions" -eq 1 ]]; then 
			
		echo " we have one session in the BIDS dir, this is good."
			
		# now we need to search for the images
		# then also find which modalities are available and set wf accordingly
			
		search_T1=($(find $search_sessions -type f | grep T1w.nii.gz));
		# search_T2=($(find $search_sessions -type f | grep T2w.nii.gz));
		# search_FLAIR=($(find $search_sessions -type f | grep FLAIR.nii.gz));
			
		if [[ $search_T1 ]]; then
				
			T1_orig=$search_T1
			echo " We found T1 WIs ${T1_orig}"
				
		else
				
			echo " no T1 WIs found in BIDS dir, exiting"
			exit 2
				
		fi

	else 
			
		echo " There is a problem with sessions in BIDS dir. "
		echo " Please double check your data structure &/or specify one session with -s if you have multiple ones. "
		exit 2
			
	fi

    if [[ "$o_flag" -eq 0 ]]; then

        output_d="${cwd}/BIDS/derivatives/output_VBG/sub-${subj}${ses_long}"

    fi

    if [[ "$m_flag" -eq 0 ]]; then

        preproc="${cwd}/BIDS/derivatives/proc_VBG/sub-${subj}${ses_long}"

    fi

elif [[ "$bids_flag" -eq 1 ]] && [[ "$s_flag" -eq 1 ]]; then
		
	# this is fine
    ses_string="${cwd}/BIDS/sub-${subj}_ses-${ses}"
	search_sessions=($(find ${ses_string} -type d | grep anat));
	num_sessions=1;
	ses_long=_ses-0${num_sessions};
		
	if [[ "$num_sessions" -eq 1 ]]; then 
			
		echo " One session " $ses " specified in BIDS dir, good."
		
		search_T1=($(find $search_sessions -type f | grep T1w.nii.gz));
		# search_T2=($(find $search_sessions -type f | grep T2w.nii.gz));
		# search_FLAIR=($(find $search_sessions -type f | grep flair.nii.gz));
		
		if [[ "$search_T1" ]]; then
			
			T1_orig=$search_T1;

            echo " We found T1 WIs ${T1_orig}"
			
		else
			
			echo " no T1 WIs found in BIDS dir, exiting "

			exit 2
			
		fi

    fi

    if [[ "$o_flag" -eq 0 ]]; then

        output_d="${cwd}/BIDS/derivatives/output_VBG/sub-${subj}${ses_long}"

    fi

    if [[ "$m_flag" -eq 0 ]]; then

        preproc="${cwd}/BIDS/derivatives/proc_VBG/sub-${subj}${ses_long}"

    fi

elif [[ "$bids_flag" -eq 0 ]] && [[ "$s_flag" -eq 0 ]]; then

	# this is fine if T1 and T2 and/or flair are set
	# find which ones are set and define wf accordingly
    num_sessions=1;
    ses_long="";
		
	if [[ "$t1_flag" ]]; then
			
		T1_orig=$t1_orig

        echo " T1 images provided as ${t1_orig} "
		
    else

        echo " No T1 WIs specified, exiting. "

		exit 2
			
	fi
		
		
elif [[ "$bids_flag" -eq 0 ]] && [[ "$s_flag" -eq 1 ]]; then
			
	echo " Wrong optional arguments, we cant have sessions without BIDS, exiting."
    
	exit 2
		
fi


######

ROIs="${output_d}/sub-${subj}${ses_long}/ROIs"
	
overlap="${output_d}/sub-${subj}${ses_long}/overlap"

#####

# make your dirs

mkdir -p ${preproc_m} >/dev/null 2>&1

mkdir -p ${output_m} >/dev/null 2>&1

mkdir -p ${preproc} >/dev/null 2>&1

mkdir -p ${output_d} >/dev/null 2>&1

mkdir -p ${ROIs} >/dev/null 2>&1

mkdir -p ${overlap} >/dev/null 2>&1


# make your log file

prep_log="${preproc}/KUL_VBG_prep_log_${d}.txt";

if [[ ! -f ${prep_log} ]] ; then

    touch ${prep_log}

else

    echo "${prep_log} already created"

fi

echo " Preproc dir is ${preproc} and output dir is ${output_d}" | tee -a ${prep_log}

echo " You are using KUL_VBG.sh version ${v}" | tee -a ${prep_log}


# deal with ncpu and itk ncpu

# itk default ncpu for antsRegistration
itk_ncpu="export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${ncpu}"
export $itk_ncpu
silent=1

# decide on BET method

if [[ -z "${BET_flag}" ]]; then

    echo
    echo " You have not specified a BET method, ANTsBET will be used by default" | tee -a ${prep_log}
    echo
    BET_m=2

else

    if [[ ${BET_m} -eq 1 ]]; then
    
        echo
        echo " You have specified HD-BET for brain extraction, please make sure it is called correctly from within KUL_VBG"
        echo " You have specified HD-BET for brain extraction, please make sure it is called correctly from within KUL_VBG" | tee -a ${prep_log}
        echo " In case of BET problems see lines 1124 - 1200 "
        echo
        # BET_m=1

    elif [[ ${BET_m} -eq 2 ]]; then

        echo
        echo " You have specified ANTs-BET for brain extraction, please make sure it is called correctly from within KUL_VBG"
        echo " You have specified ANTs-BET for brain extraction, please make sure it is called correctly from within KUL_VBG" | tee -a ${prep_log}
        echo " In case of BET problems see lines 1124 - 1200 "
        echo
        # BET_m=2

    else 

        echo
        echo " You have specified an incorrect value to the -B option, exiting... "
        echo " Correct options for the -B flag are 1 for HD-BET or 2 for ANTs-BET"
        exit 2

    fi
    
fi

# set this manually for debugging
function_path=($(which KUL_VBG.sh | rev | cut -d"/" -f2- | rev))
mrtrix_path=($(which mrmath | rev | cut -d"/" -f3- | rev))
FS_path1=($(which recon-all | rev | cut -d"/" -f3- | rev))

if [[  -z  ${function_path}  ]]; then

    echo "update function path to reflect funciton name line 514"
    exit 2

else

    echo " VBG lives in ${function_path} "
    echo " VBG lives in ${function_path} " | tee -a ${prep_log}

fi

#  the primary image is the noncontrast T1

prim=${T1_orig}

# this if loop will quit the script if the T1 is not found

if [[ -z "${T1_orig}" ]]; then

    echo
    echo " Incorrect T1 input, please check the file path and name "
    echo
    exit 2

else

    echo "Inputs are -p  ${subj}  -T1 ${T1_orig}  -lesion  ${L_mask}  -lesion_space  ${L_mask_space}"
    echo "Inputs are -p  ${subj}  -T1 ${T1_orig}  -lesion  ${L_mask}  -lesion_space  ${L_mask_space}" | tee -a ${prep_log}
    
fi


# REST OF SETTINGS ---

# Some parallelisation

if [[ "$n_flag" -eq 0 ]]; then

	ncpu=8

	echo " -n flag not set, using default 8 threads. "
    echo " -n flag not set, using default 8 threads. " | tee -a ${prep_log}

else

	echo " -n flag set, using " ${ncpu} " threads."
    echo " -n flag set, using " ${ncpu} " threads." | tee -a ${prep_log}

fi

echo "KUL_VBG @ ${d} with parent pid $$ "
echo "KUL_VBG @ ${d} with parent pid $$ " | tee -a ${prep_log}

# --- MAIN ----------------
# Start with your Vars for Part 1

# naming strings

    str_pp="${preproc}/sub-${subj}${ses_long}"

    str_op="${output_d}/sub-${subj}${ses_long}"

    str_overlap="${overlap}/sub-${subj}${ses_long}"

# Template stuff

# check which template to use based on
# is this a pediatric or adult brain and whether we use donor tissue or not

if [[ "${p_flag}" -eq 1 ]] && [[ "${t_flag}" -eq 0 ]]; then
    # ADJUST TEMPLATES FOR NKI10U IF P=1 T=0

    echo "Working with default pediatric template and priors"
    echo "Working with default pediatric template and priors" | tee -a ${prep_log}

    MNI_T1="${function_path}/atlasses/Templates_update/VBG_"

    MNI_T1_brain="${function_path}/atlasses/Templates_update/NKI10u_temp_brain.nii.gz"

    MNI_brain_mask="${function_path}/atlasses/Templates_update/NKI10u_temp_brain_mask.nii.gz"

    MNI_brain_pmask="${function_path}/atlasses/Templates_update/ped_PBEM.nii.gz"

    MNI_brain_emask="${function_path}/atlasses/Templates_update/ped_BET_mask.nii.gz"

    new_priors="${function_path}/atlasses/Templates_update/priors/NKI10U_Prior_%d.nii.gz"

elif [[ "${p_flag}" -eq 1 ]] && [[ "${t_flag}" -eq 1 ]]; then
    # ADJUST TEMPLATES FOR VBG_PED IF P=1 T=1

    echo "Working with cooked template and priors"
    echo "Working with cooked template and priors" | tee -a ${prep_log}

    MNI_T1="${function_path}/atlasses/Templates_update/VBG_T1_temp_ped.nii.gz"

    MNI_T1_brain="${function_path}/atlasses/Templates_update/VBG_T1_temp_ped_brain.nii.gz"

    MNI_brain_mask="${function_path}/atlasses/Templates_update/VBG_T1_temp_ped_brain_mask.nii.gz"

    MNI_brain_pmask="${function_path}/atlasses/Templates_update/ped_PBEM.nii.gz"

    MNI_brain_emask="${function_path}/atlasses/Templates_update/ped_BET_mask.nii.gz"

    new_priors="${function_path}/atlasses/Templates_update/priors/VBG_ped_T_Prior_%d.nii.gz"

elif [[ "${p_flag}" -eq 0 ]] && [[ "${t_flag}" -eq 1 ]]; then

    echo "Working with cooked adult template and priors"
    echo "Working with cooked adult template and priors" | tee -a ${prep_log}

    MNI_T1="${function_path}/atlasses/Templates_update/VBG_T1_temp.nii.gz"

    MNI_T1_brain="${function_path}/atlasses/Templates_update/VBG_T1_temp_brain.nii.gz"

    MNI_brain_mask="${function_path}/atlasses/Templates_update/VBG_T1_temp_brain_mask.nii.gz"

    MNI_brain_pmask="${function_path}/atlasses/Templates_update/adult_PBEM.nii.gz"

    MNI_brain_emask="${function_path}/atlasses/Templates_update/adult_BET_mask.nii.gz"

    new_priors="${function_path}/atlasses/Templates_update/priors/VBG_adult_T_Prior_%d.nii.gz"

elif [[ "${p_flag}" -eq 0 ]] && [[ "${t_flag}" -eq 0 ]]; then

    echo "Working with default adult template and priors"
    echo "Working with default adult template and priors" | tee -a ${prep_log}

    MNI_T1="${function_path}/atlasses/Templates_update/HR_T1_MNI.nii.gz"

    MNI_T1_brain="${function_path}/atlasses/Templates_update/HR_T1_MNI_brain.nii.gz"

    MNI_brain_mask="${function_path}/atlasses/Templates_update/HR_T1_MNI_brain_mask.nii.gz"

    MNI_brain_pmask="${function_path}/atlasses/Templates_update/adult_PBEM.nii.gz"

    MNI_brain_emask="${function_path}/atlasses/Templates_update/adult_BET_mask.nii.gz"

    new_priors="${function_path}/atlasses/Templates_update/priors/HRT1_Prior_%d.nii.gz"

fi



MNI2_in_T1="${str_pp}_T1_brain_inMNI2_InverseWarped.nii.gz"

MNI2_in_T1_hm="${str_pp}_T1_brain_inMNI2_InverseWarped_HistMatch.nii.gz"

MNI_r="${function_path}/atlasses/Templates_update/Rt_hemi_mask.nii.gz"

MNI_l="${function_path}/atlasses/Templates_update/Lt_hemi_mask.nii.gz"

MNI_lw="${str_pp}_MNI_L_insubjT1_inMNI1.nii.gz"

MNI_lwr="${str_pp}_MNI_L_insubjT1_inMNI1r.nii.gz"

MNI_rwr="${str_pp}_MNI_R_insubjT1_inMNI1r.nii.gz"

L_hemi_mask="${str_pp}_L_hemi_mask_bin.nii.gz"

H_hemi_mask="${str_pp}_H_hemi_mask_bin.nii.gz"

L_hemi_mask_binv="${str_pp}_L_hemi_mask_binv.nii.gz"

H_hemi_mask_binv="${str_pp}_H_hemi_mask_binv.nii.gz"

# CSF+GMC+GMB+WM  and the rest

tmp_s2T1_nCSFGMC="${str_pp}_tmp_s2T1_nCSFGMC.nii.gz"

tmp_s2T1_nCSFGMCB="${str_pp}_tmp_s2T1_nCSFGMCB.nii.gz"

tmp_s2T1_nCSFGMCBWM="${str_pp}_tmp_s2T1_nCSFGMCBWM.nii.gz"

tmp_s2T1_nCSFGMCBWMr="${str_pp}_tmp_s2T1_nCSFGMCBWMr.nii.gz"

tmp_s2T1_CSFGMCBWM="${str_pp}_tmp_s2T1_CSFGMCBWM.nii.gz"

MNI2_in_T1_scaled="${str_pp}_MNI_brain_IW_scaled.nii.gz"

tissues=("CSF" "GMC" "GMBG" "WM");

priors_str="${new_priors::${#new_priors}-9}*.nii.gz"

priors_array=($(ls ${priors_str}))

if [[ -z ${priors_array} ]]; then 

    echo " priors are not found!"
    exit 2

else


    echo "priors are ${priors_array}"

fi
    
# arrays

declare -a Atropos1_posts

declare -a Atropos2_posts

# need also to declare tpm arrays

declare -a atropos1_tpms_Lfill

declare -a atropos2_tpms_filled

declare -a atropos2_tpms_filled_GLC

declare -a atropos2_tpms_filled_GLCbinv

declare -a atropos2_tpms_punched

declare -a NP_arr_rs

declare -a NP_arr_rs_bin

declare -a NP_arr_rs_binv

declare -a NP_arr_rs_bin2

declare -a NP_arr_rs_binv2

declare -a Atropos2_posts_bin

declare -a Atropos2_posts_bin2

declare -a Atropos1_posts_bin

declare -a T1_ntiss_At2masked

declare -a nMNI2_inT1_ntiss_sc2T1MNI1

declare -a MNI2_inT1_ntiss

declare -a Atropos2_Int_finder

declare -a R_nTiss_Norm_mean

declare -a R_nTiss_Int_map_norm

declare -a Atropos1b_ntiss_map

declare -a A1_nTiss_Norm_mean

declare -a A1_nTiss_Int_scaled

declare -a A1_nTiss_Int_scaled_fill

declare -a R_nTiss_map_filled


# input variables

# lesion stuff

Lmask_o=$L_mask

L_mask_reori="${str_pp}_L_mask_reori.nii.gz"

L_mask_reori1="${str_pp}_L_mask_reori1.nii.gz"

Lmask_FS_reori="${str_pp}_Lmask_FS_reori.nii.gz"

L_mask_affMNI1="${str_pp}_L_mask_r_MNI1aff.nii.gz"

L_mask_MNI1c_bin="${str_pp}_L_mask_r_MNI1aff_bin.nii.gz"

L_mask_MNI1c_binv="${str_pp}_L_mask_r_MNI1aff_binv.nii.gz"

L_O_binv="${str_pp}_L_mask_reori_binv.nii.gz"

Lmask_bin="${str_pp}_L_mask_orig_bin.nii.gz"

Lmask_in_T1="${str_pp}_L_mask_in_T1.nii.gz"

Lmask_in_T1_bin="${str_pp}_L_mask_in_T1_bin.nii.gz"

Lmask_in_T1_binv="${str_pp}_L_mask_in_T1_binv.nii.gz"

Lmask_bin_s3="${str_pp}_Lmask_in_T1_bins3.nii.gz"

Lmask_bin_s3_flat="${str_pp}_Lmask_in_T1_bins3_flat.nii.gz"

Lmask_binv_s3_nobrain="${str_pp}_Lmask_in_T1_binvs3_nobrain.nii.gz"

Lmask_binv_s3="${str_pp}_Lmask_in_T1_binvs3.nii.gz"

brain_mask_minL="${str_pp}_antsBET_BrainMask_min_L.nii.gz"

brain_mask_minL_inMNI1="${str_pp}_brainmask_minL_inMNI1.nii.gz"

brain_mask_minL_inMNI2="${str_pp}_brainmask_minL_inMNI2.nii.gz"

brain_mask_minL_atropos2="${str_pp}_brainmask_minL_atropos2.nii.gz"

Lmask_bin_inMNI1="${str_pp}_Lmask_bin_inMNI1.nii.gz"

Lmask_binv_inMNI1="${str_pp}_Lmask_binv_inMNI1.nii.gz"

Lmask_bin_inMNI1_s3="${str_pp}_Lmask_bin_s3_inMNI1.nii.gz"

Lmask_bin_inMNI1_dilx2="${str_pp}_Lmask_bin_inMNI1_dilmx2.nii.gz"

Lmask_binv_inMNI1_dilx2="${str_pp}_Lmask_binv_inMNI1_dilmx2.nii.gz"

Lmask_binv_inMNI1_s3="${str_pp}_Lmask_binv_s3_inMNI1.nii.gz"

Lmask_bin_inMNI2_s3="${str_pp}_Lmask_bin_s3_inMNI2.nii.gz"

Lmask_binv_inMNI2_s3="${str_pp}_Lmask_binv_s3_inMNI2.nii.gz"

Lmask_bin_inMNI2="${str_pp}_Lmask_bin_inMNI2.nii.gz"

fbrain_mask_minL_inMNI1="${str_pp}_fbrainmask_minL_inMNI1.nii.gz"

L_fill_T1="${str_pp}_T1_Lfill_inMNI2.nii.gz"

nat_T1_filled1="${str_pp}_T1inMNI2_fill1.nii.gz"

stitched_T1_temp="${str_pp}_stitched_T1_brain_temp.nii.gz"

stitched_T1_nat="${str_pp}_stitched_T1_brain_nat.nii.gz"

stitched_T1_innat="${str_pp}_stitched_T1_brain_bk2nat.nii.gz"

# stitched_T1_temp_innat="${str_pp}_stitched_T1_brain_temp_bk2nat.nii.gz"

# stitched_T1="${str_pp}_stitched_T1_brain.nii.gz"

T1_bk2nat1_str="${str_pp}_T1_brain_bk2anat1_"

Temp_L_hemi="${str_pp}_Temp_L_hemi_filler.nii.gz"

Temp_L_fill_T1="${str_pp}_Temp_L_fill_T1.nii.gz"

Temp_T1_bilfilled1="${str_pp}_T1_brain_Temp_bil_Lmask_filled1.nii.gz"

Temp_bil_Lmask_fill1="${str_pp}_Temp_bil_Lmask_fill1.nii.gz"

Temp_T1_filled1="${str_pp}_Temp_T1inMNI2_filled1.nii.gz"

T1_filled_bk2nat1="${str_pp}_T1_brain_bk2anat1_InverseWarped.nii.gz"

filled_segm_im="${str_pp}_atropos1_Segmentation_2nat.nii.gz"

real_segm_im="${str_pp}_atropos2_Segmentation_2nat.nii.gz"

Lfill_segm_im="${str_pp}_Lfill_segmentation_im.nii.gz"

atropos2_segm_im_filled="${str_pp}_filled_atropos2_segmentation_im.nii.gz"

atropos2_segm_im_filled_nat="${str_pp}_filled_atropos2_segmentation_im.nii.gz"

lesion_left_overlap="${str_overlap}_L_lt_overlap.nii.gz"

lesion_right_overlap="${str_overlap}_L_rt_overlap.nii.gz"

smoothed_binLmask15="${str_pp}_smoothedLmaskbin15.nii.gz"

smoothed_binvLmask15="${str_pp}_smoothedLmaskbinv15.nii.gz"

# last lesion related vars (hopefully)

L_mask_reori_scaled="${str_pp}_L_mask_reori_scaled99.nii.gz"

bmc_minL_conn="${str_pp}_brain_mask_cleaned_minL_conn.nii.gz"

bmc_minL_true="${str_pp}_brain_mask_cleaned_minL.nii.gz"

L_mask_reori_ero1="${str_pp}_L_mask_reori_ero1.nii.gz"

bmc_minL_ero1="${str_pp}_brain_mask_cleaned_minL_ero1.nii.gz"

L_mask_reori_ero2="${str_pp}_L_mask_reori_ero2.nii.gz"

bmc_minL_ero2="${str_pp}_brain_mask_cleaned_minL_ero2.nii.gz"

# img vars for part 1 and 2

T1_reori_mat="${str_pp}_T1_reori2std_matrix.mat"

T1_reori_mat_inv="${str_pp}_T1_reori2std_matrix_inv.mat"

# T1_pp1="${str_pp}_T1_bfc.nii.gz"

T1_pp1="${str_pp}_T1_dn_thr.nii.gz"

T1_brain="${str_pp}_antsBET_BrainExtractionBrain.nii.gz"

brain_mask="${str_pp}_antsBET_BrainExtractionMask.nii.gz"

T1_inMNI_aff_str="${str_pp}_T1_inMNI_aff"

T1_inMNI_aff="${str_pp}_T1_inMNI_aff_Warped.nii.gz"

KULBETp="${str_pp}_atropos4BET"

rough_mask="${str_pp}_rough_mask.nii.gz"

rough_mask_minL="${str_pp}_rough_mask_minL.nii.gz"

clean_mask_nat_binv="${str_pp}_clean_brain_mask_nat_binv.nii.gz"

T1_brain_clean="${str_pp}_Brain_clean.nii.gz"

MNI_bm_BET_innat="${str_pp}_MNI_BM_inNat.nii.gz"

clean_mask_nat="${str_pp}_Brain_clean_mask.nii.gz"

clean_BM_mgz="${str_pp}_Brain_clean_mask.mgz"

hdbet_str="${str_pp}_Brain_clean"

BET_mask_s2="${str_pp}_antsBET_Mask_s2.nii.gz"

BET_mask_binvs2="${str_pp}_antsBET_Mask_binv_s2.nii.gz"

T1_skull="${str_pp}_T1_skull.nii.gz"

T1_brMNI1_str="${str_pp}_T1_brain_inMNI1_"

T1_brain_inMNI1="${str_pp}_T1_brain_inMNI1_Warped.nii.gz"

# T1_noise_inMNI1="${str_pp}_T1_noise_inMNI1.nii.gz"

# fT1_noise_inMNI1="${str_pp}_fT1_noise_inMNI1.nii.gz"

# T1_noise_H_hemi="${str_pp}_T1_noise_Hhemi_inMNI1.nii.gz"

# stitched_noise_MNI1="${str_pp}_T1_stitched_noise_inMNI1.nii.gz"

# stitched_noise_nat="${str_pp}_T1_stitched_noise_nat.nii.gz"

T1_brMNI2_str="${str_pp}_T1_brain_inMNI2_"

T1_brain_inMNI2="${str_pp}_T1_brain_inMNI2_Warped.nii.gz"

fT1brain_inMNI1="${str_pp}_fT1_brain_inMNI1_Warped.nii.gz"

fT1_brMNI2_str="${str_pp}_fT1_brain_inMNI2_"

fT1brain_inMNI2="${str_pp}_fT1_brain_inMNI2_Warped.nii.gz"

brain_mask_inMNI1="${str_pp}_brain_mask_inMNI1.nii.gz"

# MNI_brain_mask_in_nat="${str_pp}_MNI_brain_mask_in_nat.nii.gz"

T1_sti2fil_str="${str_pp}_stitchT12filled_brain_"

T1fill2MNI1minL_str="${str_pp}_filledT12MNI1_brain_"

T1_sti2fill_brain="${str_pp}_stitchT12filled_brain_Warped.nii.gz"

T1_fin_Lfill_1="${str_pp}_T1_finL_fill_1.nii.gz"

T1_fin_Lfill_2="${str_pp}_T1_finL_fill_2.nii.gz"

T1_fin_filled="${str_pp}_T1_finL_filled.nii.gz"

Lmask_binv_s3_n_ori="${str_pp}_Lmask_binv_s3_orig_ori.nii.gz"

T1_fin_Lfill_n_ori="${str_pp}_T1_final_Lesion_fill.nii.gz"

T1_nat_filled_out_1="${str_pp}_T1_stdOri_filled_1.nii.gz"

T1_nat_fout_wskull_1="${str_pp}_T1_stdOri_filld_wskull_1.nii.gz"

T1_nat_filled_out_2="${str_op}_T1_stdOri_filled.nii.gz"

T1_nat_fout_wskull_2="${str_op}_T1_stdOri_filld_wskull.nii.gz"

# vars for final output in input space

T1_4_FS="${str_op}_T1_nat_filled.nii.gz"

T1_4_parc="${str_op}_T1_nat_4parc.mgz"

T1_Brain_4_FS="${str_op}_T1_nat_filled_brain.nii.gz"

T1_BM_4_FS="${str_op}_T1_nat_filled_mask.nii.gz"

# img vars for part 2

T1_H_hemi="${str_pp}_T1_H_hemi.nii.gz"

fT1_H_hemi="${str_pp}_fT1_H_hemi.nii.gz"

# img vars for make ims loops

MNI2_in_T1_linsc_norm="${str_pp}_MNI2_inT1_linsc_norm.nii.gz"

atropos1_brain_norm="${str_pp}_Atropos1_brain_norm.nii.gz"

T1b_inMNI1_pN_sc2st2f="${str_pp}_T1b_inMNI1_pN_sc2_st2fill.nii.gz"

T1b_inMNI1_punched="${str_pp}_T1brain_inMNI1_punched.nii.gz"

T1b_inMNI1_p_norm="${str_pp}_T1brain_inMNI1_punched_norm.nii.gz"

#

# workflow markers for processing control

search_wf_mark1=($(find ${preproc} -type f | grep "${brain_mask_inMNI1}"));

srch_preprocp1=($(find ${preproc} -type f | grep "${T1_pp1}"));

srch_antsBET=($(find ${preproc} -type f | grep "${T1_brain_clean}"));

T1brain2MNI1=($(find ${preproc} -type f | grep "${T1_brain_inMNI1}"));

T1brain2MNI2=($(find ${preproc} -type f | grep "${T1_brain_inMNI2}"));

fT1_brain_2MNI2=($(find ${preproc} -type f | grep "${fT1brain_inMNI2}"));

sch_brnmsk_minL=($(find ${preproc} -type f | grep "${brain_mask_minL}"));

search_wf_mark2=($(find ${preproc} -type f | grep "${fbrain_mask_minL_inMNI1}"));

srch_Lmask_pt2=($(find ${preproc} -type f | grep "*_atropos1_Segmentation.nii.gz")); # same as Atropos1 marker for now

stitch2fill_mark=($(find ${preproc} -type f | grep "${T1_sti2fil_str}Warped.nii.gz"));

fill2MNI1mniL_mark=($(find ${preproc} -type f | grep "${T1fill2MNI1minL_str}Warped.nii.gz"));

Atropos1_wf_mark=($(find ${preproc} -type f | grep "_atropos1_Segmentation.nii.gz"));

srch_bk2anat1_mark=($(find ${preproc} -type f | grep "${T1_filled_bk2nat1}"));

Atropos2_wf_mark=($(find ${preproc} -type f | grep "_atropos2_Segmentation.nii.gz"));

srch_postAtropos2=($(find ${preproc} -type f | grep "_atropos2_SegmentationPosterior2_clean_bk2nat1.nii.gz"));

srch_make_images=($(find ${output_d} -type f | grep "${T1_BM_4_FS}")); # search make images will run at the end if the BM for FS is not found

# Misc subfuctions

# execute function (maybe add if loop for if silent=0)

function task_exec {

    echo "  " | tee -a ${prep_log} 
    
    echo ${task_in} | tee -a ${prep_log} 

    echo " Started @ $(date "+%Y-%m-%d_%H-%M-%S")" | tee -a ${prep_log} 

    eval ${task_in} | tee -a ${prep_log} 2>&1 &

    echo " pid = $! basicPID = $$ " | tee -a ${prep_log}

    wait ${pid}

    sleep 5

    if [ $? -eq 0 ]; then
        echo Success | tee -a ${prep_log}
    else
        echo Fail | tee -a ${prep_log}

        exit 1
    fi

    echo " Finished @  $(date "+%Y-%m-%d_%H-%M-%S")" | tee -a ${prep_log} 

    echo "  " | tee -a ${prep_log} 

    unset task_in

}

# functions for basic antsRegSyN calls

# not using SyNQuick anymore
# default Affine antsRegSyNQuick call
# function KUL_antsRegSyNQ_Def {

#     task_in="antsRegistrationSyNQuick.sh -d 3 -f ${fix_im} -m ${mov_im} -o ${output} -n ${ncpu} -j 1 -t ${transform} ${mask}"

#     task_exec

# }

# default Affine antsRegSyN call
function KUL_antsRegSyN_Def {

    task_in="antsRegistrationSyN.sh -d 3 -f ${fix_im} -m ${mov_im} -o ${output} -n ${ncpu} -j 1 -t ${transform} ${mask}"

    task_exec

}

# functions for ANTsBET

# adding new ANTsBET workflow
# actually, we could use hd-bet cpu version if it is installed also
# make a little if loop testing if hd-bet is alive

function KUL_antsBETp {

    # task_in="fslreorient2std ${Lmask_o} ${L_mask_reori}"

    # task_exec

    # BET is done after an initial affine transform to template space

    task_in="antsRegistrationSyN.sh -d 3 -f ${MNI_T1} -m ${prim_in} -o ${output}_aff_2_temp_ -t a"
    
    task_exec

    task_in="antsApplyTransforms -d 3 -i ${L_mask_reori} -o ${L_mask_affMNI1} -r ${MNI_T1} -t [${output}_aff_2_temp_0GenericAffine.mat,0] \
    && fslmaths ${L_mask_affMNI1} -mas ${MNI_brain_mask} -bin -save ${L_mask_MNI1c_bin} -binv ${L_mask_MNI1c_binv}"

    task_exec

    # this approach ensures minimal failures in either case
    # if HD-BET excludes too much of a brain or if ANTs includes too much

    if [[ ${BET_m} -eq 1 ]]; then

        echo "HD-BET is selected, will use this for brain extraction" | tee -a ${prep_log}

        echo "sourcing ptc conda virtual env, if yours is named differently please edit lines 822 823 " | tee -a ${prep_log}

        # task_in="source /anaconda3/bin/activate ptc && hd-bet -i ${prim_in} -o ${output} -tta 0 -mode fast -s 1 -device cpu"

        # task_exec

        task_in="hd-bet -i ${output}_aff_2_temp_Warped.nii.gz -o ${output}_i -tta 0 -mode fast -s 1 -device cpu"

        task_exec

        task_in="mrcalc -force -nthreads ${ncpu} ${output}_i_mask.nii.gz ${L_mask_MNI1c_binv} -mul ${L_mask_MNI1c_bin} -add ${output}_brain_mask_c_h_MNI1aff.nii.gz \
        && ImageMath 3 ${output}_brain_mask_c_hf_MNI1aff.nii.gz FillHoles ${output}_brain_mask_c_h_MNI1aff.nii.gz \
        && mrcalc -force -nthreads ${ncpu} ${output}_brain_mask_c_hf_MNI1aff.nii.gz ${output}_aff_2_temp_Warped.nii.gz -mult ${output}_brain_c_MNI1aff.nii.gz"

        task_exec

        task_in="antsApplyTransforms -d 3 -i ${output}_brain_c_MNI1aff.nii.gz -o ${T1_brain_clean} -r ${prim_in} -t [${output}_aff_2_temp_0GenericAffine.mat,1] \
        && antsApplyTransforms -d 3 -i ${output}_brain_mask_c_hf_MNI1aff.nii.gz -o ${output}_brain_mask_clean_innat_lin.nii.gz -r ${prim_in} -t [${output}_aff_2_temp_0GenericAffine.mat,1] -n MultiLabel\
        && fslmaths ${output}_brain_mask_clean_innat_lin.nii.gz -bin ${clean_mask_nat}"
        # && ImageMath 3 ${output}_brain_mask_clean_innat_lin_FH.nii.gz FillHoles ${output}_brain_mask_clean_innat_lin.nii.gz \
        

        task_exec

    elif [[ ${BET_m} -eq 2 ]]; then

        echo "ANTsBET is selected, will use this for brain extraction" | tee -a ${prep_log}

        # if you want to use the modified ANTs based BET approach and not HD-BET
        # just comment out the if loop and hd-bet condition (be sure to get the if, else and fi lines)
        
        task_in="antsBrainExtraction.sh -d 3 -a ${output}_aff_2_temp_Warped.nii.gz -e ${MNI_T1} -m ${MNI_brain_pmask} -f ${MNI_brain_emask} -u 1 -k 1 -q 1 -o ${output}_"

        task_exec

        # Bring results back to native space

        task_in="fslmaths ${output}_BrainExtractionMask.nii.gz -mul ${MNI_brain_pmask} -save ${output}_brain_mask_c_MNI1aff.nii.gz -restart \
        ${output}_BrainExtractionBrain.nii.gz -mul ${output}_brain_mask_c_MNI1aff.nii.gz ${output}_BrainExtractionBrain_c.nii.gz \
        && antsApplyTransforms -d 3 -i ${output}_BrainExtractionBrain_c.nii.gz -o ${T1_brain_clean} -r ${prim_in} -t [${output}_aff_2_temp_0GenericAffine.mat,1]"

        task_exec

        task_in="fslmaths ${T1_brain_clean} -bin ${clean_mask_nat}"

        task_exec

    else

        echo "we have a problem"

    fi

    # exit 2

}


# Dealing with the lesion mask part 1

# function KUL_Lmask_part1 {

#     # since we only operate in 1 space (unimodal) this if condition is useless and deprecated
#     # substituting with E_flag coniditional arguments

#     if [[ "${E_flag}" -eq 0 ]]; then

#         # echo " Lesion mask is already in T1 space " | tee -a ${prep_log}

#         echo " Intra-axial lesion running VBG Lmask_pt1 workflow and subsequent steps" | tee -a ${prep_log}

#         # start by smoothing and thring the mask

#         task_in="fslmaths ${L_mask_reori} -s 2 -thr 0.2 -bin -save ${Lmask_bin} -binv ${Lmask_in_T1_binv}"

#         task_exec

#         echo " Copying Lmask_bin_s2 file to Lmask_in_T1_bin " | tee -a ${prep_log}

#         cp ${Lmask_bin} ${Lmask_in_T1_bin}

#         # subtract lesion from brain mask

#         task_in="fslmaths ${clean_mask_nat} -mas ${Lmask_in_T1_binv} -mas ${clean_mask_nat} ${brain_mask_minL}"

#         task_exec

#     else

#         echo " Extra-axial lesion running simplified VBG Lmask_pt1 workflow, FS and subsequent steps" | tee -a ${prep_log}

#         task_in="fslmaths ${L_mask_reori} -binv ${L_O_binv}"

#         task_exec

    
#     fi

# }

#  determine lesion laterality and proceed accordingly
#  define all vars for this function

function KUL_Lmask_part2 {
    
    #############################################################
    # creating edited lesion masks
    # should add proc control to this section

    # now do the unflipped brain_mask_minL & L_mask

    mask_in="${Lmask_bin_inMNI1}"

    mask_out="${Lmask_bin_inMNI2}"

    ref="${MNI_T1_brain}"

    task_in="antsApplyTransforms -d 3 -i ${mask_in} -o ${mask_out} -r ${ref} -t ${T1_brMNI2_str}1Warp.nii.gz -t [${T1_brMNI2_str}0GenericAffine.mat,0] -n MultiLabel"

    task_exec

    unset mask_in mask_out

    task_in="fslmaths ${Lmask_bin_inMNI2} -bin -save ${Lmask_bin_inMNI2} -binv -mas ${brain_mask_inMNI1} ${brain_mask_minL_inMNI2}"

    task_exec

    # here we are really making the Lmask bigger
    # using -dilM x2 and -s 2 with a -thr 0.2 to avoid very low value voxels
    # actually will not be using dilx2 probably

    task_in="fslmaths ${Lmask_bin_inMNI1} -binv -mas ${brain_mask_inMNI1} -save ${Lmask_binv_inMNI1} -restart ${Lmask_bin_inMNI1} -dilM -dilM -save ${Lmask_bin_inMNI1_dilx2} \
    -s 2 -thr 0.2 -mas ${brain_mask_inMNI1} ${Lmask_bin_inMNI1_s3} && fslmaths ${brain_mask_inMNI1} -sub ${Lmask_bin_inMNI1_s3} -mas ${brain_mask_inMNI1} \
    ${Lmask_binv_inMNI1_s3} && fslmaths ${Lmask_bin_inMNI1_dilx2} -binv ${Lmask_binv_inMNI1_dilx2}"

    task_exec

    # task_in="fslmaths ${Lmask_bin_inMNI2} -dilM -dilM -save ${Lmask_bin_inMNI2_dilx2} -s 2 -thr 0.2 -mas ${brain_mask_inMNI1} ${Lmask_bin_inMNI2_s3} && fslmaths \
    # ${brain_mask_inMNI1} -sub ${Lmask_bin_inMNI2_s3} -mas ${brain_mask_inMNI1} ${Lmask_binv_inMNI2_s3} && fslmaths ${Lmask_bin_inMNI1_dilx2} -binv ${Lmask_binv_inMNI1_dilx2}"

    task_in="fslmaths ${Lmask_bin_inMNI2} -dilM -dilM -s 2 -thr 0.2 -mas ${brain_mask_inMNI1} ${Lmask_bin_inMNI2_s3} && fslmaths \
    ${brain_mask_inMNI1} -sub ${Lmask_bin_inMNI2_s3} -mas ${brain_mask_inMNI1} ${Lmask_binv_inMNI2_s3}"

    task_exec

    # task_in="fslmaths ${Lmask_in_T1_bin} -dilM -dilM -save ${Lmask_bin_dilx2} -s 2 -thr 0.2 -mas ${clean_mask_nat} ${Lmask_bin_s3} && fslmaths \
    # ${clean_mask_nat} -sub ${Lmask_bin_s3} -mas ${clean_mask_nat} ${Lmask_binv_s3} && fslmaths ${Lmask_bin_dilx2} -binv ${Lmask_binv_dilx2}"

    task_in="fslmaths ${Lmask_in_T1_bin} -dilM -dilM -s 2 -thr 0.2 -mas ${clean_mask_nat} ${Lmask_bin_s3} && fslmaths \
    ${clean_mask_nat} -sub ${Lmask_bin_s3} -mas ${clean_mask_nat} ${Lmask_binv_s3}"

    task_exec

    ######################################################################

    echo " Now running lesion magic part 2 "

    # determine lesion laterality
    # this is all happening in MNI space (between the first and second warped images)
    # apply warps to MNI_rl to match patient better
    # generate L_hemi+L_mask & H_hemi_minL_mask for unilateral lesions
    # use those to generate stitched image

    task_in="antsApplyTransforms -d 3 -i ${MNI_l} -o ${MNI_lw} -r ${T1_brain_inMNI1} -t [${T1_brMNI2_str}0GenericAffine.mat,1] -t ${T1_brMNI2_str}1InverseWarp.nii.gz -n MultiLabel \
    && fslmaths ${MNI_lw} -thr 0.1 -bin -mas ${brain_mask_inMNI1} -save ${MNI_lwr} -binv -mas ${brain_mask_inMNI1} ${MNI_rwr}"

    task_exec

    task_in="fslmaths ${MNI_lwr} -mas ${Lmask_bin_inMNI1} ${lesion_left_overlap}"

    task_exec

    task_in="fslmaths ${MNI_rwr} -mas ${Lmask_bin_inMNI1} ${lesion_right_overlap}"

    task_exec

    Lmask_tot_v=$(mrstats -force -nthreads ${ncpu} ${Lmask_bin_inMNI1} -output count -quiet -ignorezero)

    overlap_left=$(mrstats -force -nthreads ${ncpu} ${lesion_left_overlap} -output count -quiet -ignorezero)

    overlap_right=$(mrstats -force -nthreads ${ncpu} ${lesion_right_overlap} -output count -quiet -ignorezero)

    L_ovLt_2_total=$(echo ${overlap_left}*100/${Lmask_tot_v} | bc)

    L_ovRt_2_total=$(echo ${overlap_right}*100/${Lmask_tot_v} | bc)

    echo " total lesion vox count ${Lmask_tot_v}" | tee -a ${prep_log}

    echo " ov_left is ${overlap_left}" | tee -a ${prep_log}

    echo " ov right is ${overlap_right}" | tee -a ${prep_log}

    echo " ov_Lt to total is ${L_ovLt_2_total}" | tee -a ${prep_log}

    echo " ov_Rt to total is ${L_ovRt_2_total}" | tee -a ${prep_log}

    # we set a hard-coded threshold of 65, if unilat. then native heatlhy hemi is used
    # if bilateral by more than 35, template brain is used
    # # this needs to be modified, also need to include simple lesion per hemisphere overlap with percent to total hemi volume
    # this will enable us to use template or simple filling and derive mean values per tissue class form another source (as we are currently using the original images).
    # AR 09/02/2020
    # here we also need to make unilateral L masks, masked by hemi mask to overcome midline issue
    
    if [[ "${L_ovLt_2_total}" -gt 65 ]]; then

        # instead of simply copying
        # here we add the lesion patch to the L_hemi
        # task_in="cp ${MNI_l} ${L_hemi_mask}"

        task_in="fslmaths ${MNI_lwr} -add ${Lmask_bin_inMNI1} -bin -mas ${brain_mask_inMNI1} -save ${L_hemi_mask} -binv -mas ${brain_mask_inMNI1} ${L_hemi_mask_binv}"

        task_exec

        task_in="fslmaths ${MNI_rwr} -mas ${Lmask_binv_inMNI1} -bin -mas ${brain_mask_inMNI1} -save ${H_hemi_mask} -binv -mas ${brain_mask_inMNI1} ${H_hemi_mask_binv}"

        task_exec

        echo ${L_ovLt_2_total} | tee -a ${prep_log}

        echo " This patient has a left sided or predominantly left sided lesion " | tee -a ${prep_log}

        echo "${L_hemi_mask}" | tee -a ${prep_log}

        echo "${H_hemi_mask}" | tee -a ${prep_log}

        # for debugging will set this to -lt 10
        # should change it back to -gt 65 ( on the off chance you will need the bilateral condition, that is also been tested now)
        
    elif [[ "${L_ovRt_2_total}" -gt 65 ]]; then

        # task_in="cp ${MNI_r} ${L_hemi_mask}"

        task_in="fslmaths ${MNI_rwr} -add ${Lmask_bin_inMNI1} -bin -save ${L_hemi_mask} -binv ${L_hemi_mask_binv}"

        task_exec

        # task_in="cp ${MNI_l} ${H_hemi_mask}"

        task_in="fslmaths ${MNI_lwr} -add ${Lmask_bin_inMNI1} -bin -save ${H_hemi_mask} -binv ${H_hemi_mask_binv}"

        task_exec

        echo ${L_ovRt_2_total} | tee -a ${prep_log}

        echo " This patient has a right sided or predominantly right sided lesion " | tee -a ${prep_log}

        echo "${L_hemi_mask}" | tee -a ${prep_log}

        echo "${H_hemi_mask}" | tee -a ${prep_log}
        
    else 

        bilateral=1

        echo " This is a bilateral lesion with ${L_ovLt_2_total} left side and ${L_ovRt_2_total} right side, using Template T1 to derive lesion fill patch. "  | tee -a ${prep_log}

        echo " note Atropos1 will use the filled images instead of the stitched ones "  | tee -a ${prep_log}

    fi

    ################

    # for loop apply warp to each tissue tpm
    # first we scale the images then
    # histogram match within masks

    # mrstats to get the medians to match with mrcalc

    med_tmp=$(mrstats -quiet -nthreads ${ncpu} -output median -ignorezero -mask ${MNI_brain_mask} ${MNI2_in_T1})

    med_nat=$(mrstats -quiet -nthreads ${ncpu} -output median -ignorezero -mask ${brain_mask_minL_inMNI2} ${T1_brain_inMNI1})

    task_in="mrcalc -force -nthreads ${ncpu} ${med_nat} ${med_tmp} -divide ${MNI2_in_T1} -mult - | mrhistmatch linear - ${T1_brain_inMNI1} ${MNI2_in_T1_scaled} -force \
    -nthreads ${ncpu} -mask_target ${brain_mask_minL_inMNI2} -mask_input ${MNI_brain_mask}"

    task_exec

    # first we get the mean and normalized MNI2_inT1 after the scaling process done above

    MNI2_inT1_sc_mean=$(fslstats ${MNI2_in_T1_scaled} -M)

    task_in="fslmaths ${MNI2_in_T1_scaled} -div ${MNI2_inT1_sc_mean} ${MNI2_in_T1_linsc_norm}"

    task_exec

    # same for the T1brain in MNI1 image

    # before we get any values for the T1 we are trying to fill
    # we have to mask out the lesion

    # we get normalized target and source T1 images by dividing each by its mean

    task_in="fslmaths ${T1_brain_inMNI1} -mas ${Lmask_binv_inMNI1_dilx2} ${T1b_inMNI1_punched}"

    task_exec

    T1b_inMNI1p_mean=$(fslstats ${T1b_inMNI1_punched} -M)

    task_in="fslmaths ${T1b_inMNI1_punched} -div ${T1b_inMNI1p_mean} ${T1b_inMNI1_p_norm}"

    task_exec

    # for each tissue type
    # attempting to minimize intensity difference between donor and recipient images

    for ts in ${!tissues[@]}; do

        NP_arr_rs[$ts]="${str_pp}_atroposP_${tissues[$ts]}_rs.nii.gz"

        NP_arr_rs_bin[$ts]="${str_pp}_atroposP_${tissues[$ts]}_rs_bin.nii.gz"

        NP_arr_rs_binv[$ts]="${str_pp}_atroposP_${tissues[$ts]}_rs_binv.nii.gz"

        Atropos2_posts_bin[$ts]="${str_pp}_atropos2_${tissues[$ts]}_bin.nii.gz"

        MNI2_inT1_ntiss[$ts]="${str_pp}_MNItmp_${tissues[$ts]}_IM.nii.gz"

        T1_ntiss_At2masked[$ts]="${str_pp}_T1inMNI2_${tissues[$ts]}_IM.nii.gz"

        nMNI2_inT1_ntiss_sc2T1MNI1[$ts]="${str_pp}_nMNI2_inT1_linsc_norm_n${tissues[$ts]}.nii.gz"

        # warp the tissues to T1_brain_inMNI1 (first deformation)

        task_in="antsApplyTransforms -d 3 -i ${priors_array[$ts]} -o ${NP_arr_rs[$ts]} -r ${T1_brain_inMNI1} \
        -t [${T1_brMNI2_str}0GenericAffine.mat,1] -t ${T1_brMNI2_str}1InverseWarp.nii.gz"

        task_exec

        task_in="mrcalc -force -nthreads ${ncpu} ${NP_arr_rs[$ts]} 0.1 -ge ${str_pp}_atropos_${tissues[$ts]}_rs_thr.nii.gz"

        task_exec

        task_in="mrcalc -force -nthreads ${ncpu} ${str_pp}_atropos_${tissues[$ts]}_rs_thr.nii.gz 0.2 -ge 1 $((ts+1)) -replace ${NP_arr_rs_bin[$ts]} && mrcalc -force -nthreads ${ncpu} ${NP_arr_rs_bin[$ts]} -1 -mult 0 -ge ${NP_arr_rs_binv[$ts]}"

        task_exec

    done

    # combine tpms to 1 3D segmentation map

    task_in="mrcalc -force -nthreads ${ncpu} ${NP_arr_rs[3]} 0.1 -ge 1 4 -replace ${NP_arr_rs_bin[3]} && mrcalc -force -nthreads ${ncpu} ${NP_arr_rs_bin[3]} -1 -mult 0 -gt ${NP_arr_rs_binv[3]}"
    
    task_exec

    task_in="mrcalc -force -nthreads ${ncpu} ${NP_arr_rs_bin[3]} ${NP_arr_rs_binv[1]} -mult ${NP_arr_rs_binv[2]} -mult ${NP_arr_rs_binv[0]} -mult ${str_pp}_atropos_priors_WM_prepped.nii.gz"

    task_exec
    
    # task_in="fslmaths ${NP_arr_rs[2]} -thr 0.2 -save ${str_pp}_atropos_${tissues[2]}_rs_thr.nii.gz -thr 0.2 -bin -mul 3 -save ${NP_arr_rs_bin[2]} -binv -save ${NP_arr_rs_binv[2]} \
    # -restart ${NP_arr_rs[1]} -thr 0.2 -save ${str_pp}_atropos_${tissues[1]}_rs_thr.nii.gz -thr 0.2 -bin -mul 2 -save ${NP_arr_rs_bin[1]} -binv -save ${NP_arr_rs_binv[1]} -restart  \
    # ${NP_arr_rs[0]} -thr 0.2 -save ${str_pp}_atropos_${tissues[0]}_rs_thr.nii.gz -thr 0.2 -bin -mul 1 -save ${NP_arr_rs_bin[0]} -binv -save ${NP_arr_rs_binv[0]} -restart \
    # ${NP_arr_rs[3]} -thr 0.05 -save ${str_pp}_atropos_${tissues[3]}_rs_thr.nii.gz -thr 0.05 -bin -mul 4 -save ${NP_arr_rs_bin[3]} -binv -save ${NP_arr_rs_binv[3]} -restart \
    # ${NP_arr_rs_bin[3]} -mul ${NP_arr_rs_bin[3]} -mul ${NP_arr_rs_binv[1]} -mul ${NP_arr_rs_binv[2]} -mul ${NP_arr_rs_binv[0]} \
    # ${str_pp}_atropos_priors_WM_prepped.nii.gz"

    # task_exec

    # add them all up to make the segmentatiom image

    task_in="fslmaths ${NP_arr_rs_bin[1]} -mul ${NP_arr_rs_binv[2]} -add ${NP_arr_rs_bin[2]} -save ${str_pp}_atropos_priors_GMC+BG_prepped.nii.gz \
    -mul ${NP_arr_rs_binv[0]} -add ${NP_arr_rs_bin[0]} ${str_pp}_atropos_priors_GMC+BG+CSF_prepped.nii.gz && ImageMath 3 \
    ${str_pp}_atropos_priors_GMC+BG+CSF+WM_p.nii.gz addtozero ${str_pp}_atropos_priors_GMC+BG+CSF_prepped.nii.gz ${NP_arr_rs_bin[3]} \
    ${str_pp}_atropos_priors_GMC+BG+CSF+WM_p.nii.gz && ImageMath 3 ${str_pp}_atropos_priors_GMC+BG+CSF+WM_p2.nii.gz \
    ReplaceVoxelValue ${str_pp}_atropos_priors_GMC+BG+CSF+WM_p.nii.gz 5 10 4 && ImageMath 3 ${str_pp}_atropos_priors_GMC+BG+CSF+WM_ready.nii.gz \
    FillHoles ${str_pp}_atropos_priors_GMC+BG+CSF+WM_p2.nii.gz"

    task_exec

    for ts in ${!tissues[@]}; do

        NP_arr_rs_bin2[$ts]="${str_pp}_atroposP_${tissues[$ts]}_rs_bin2.nii.gz"

        NP_arr_rs_binv2[$ts]="${str_pp}_atroposP_${tissues[$ts]}_rs_binv2.nii.gz"

        Atropos2_posts_bin2[$ts]="${str_pp}_atropos2_${tissues[$ts]}_bin2.nii.gz"
       
        # create the tissue masks and punch the lesion out of them

        task_in="fslmaths ${str_pp}_atropos_priors_GMC+BG+CSF+WM_ready.nii.gz -thr $((ts+1)) -uthr $((ts+1)) -bin -mul ${NP_arr_rs[$ts]} -mas ${MNI_brain_mask} -bin -save ${NP_arr_rs_bin2[$ts]} -binv \
        ${NP_arr_rs_binv2[$ts]} && mrthreshold -force -quiet -nthreads ${ncpu} -percentile 99.5 ${Atropos2_posts[$ts]} - | mrcalc - ${brain_mask_minL_inMNI1} -mult 0.5 -ge ${Atropos2_posts_bin2[$ts]} -force"

        task_exec

    done

    task_in="fslmaths ${str_pp}_atropos_priors_GMC+BG+CSF+WM_ready.nii.gz -thr 4 -uthr 4 -bin -mul ${NP_arr_rs[3]} -mas ${MNI_brain_mask} -bin -save ${str_pp}_atroposP_WM_rs_bin2.nii.gz -binv \
    ${str_pp}_atroposP_WM_rs_binv2.nii.gz && mrthreshold -force -quiet -nthreads ${ncpu} -percentile 95 ${Atropos2_posts[3]} - | mrcalc - ${brain_mask_minL_inMNI1} -mult 0.2 -ge \
    ${str_pp}_atropos2_WM_bin2.nii.gz -force"

    task_exec

    for gs in ${!tissues[@]}; do

        # get tissue intensity maps from template

        task_in="fslmaths ${MNI2_in_T1_linsc_norm} -mas ${NP_arr_rs_bin2[$gs]} ${MNI2_inT1_ntiss[$gs]} && fslmaths ${T1b_inMNI1_p_norm} -mas ${Atropos2_posts_bin2[$gs]} \
        ${T1_ntiss_At2masked[$gs]}"

        task_exec

        T1_ntiss_At2m_mean=$(fslstats ${T1_ntiss_At2masked[$gs]} -M)

        MNI2_inT1_ntiss_mean=$(fslstats ${MNI2_inT1_ntiss[$gs]} -M)

        task_in="fslmaths ${MNI2_inT1_ntiss[$gs]} -div ${MNI2_inT1_ntiss_mean} -mul ${T1_ntiss_At2m_mean} ${nMNI2_inT1_ntiss_sc2T1MNI1[$gs]}"

        task_exec
        
    done

    # Sum up the tissues while masking in and out to minimize overlaps and holes
    # here we use ImageMath addtozero instead, to preserve a smooth interface between the tissues
    # the fslmaths step works well, but results a rather cartoon looking image
    # the order of images input to ImageMath addtozero makes a difference
    # we start with GMC as that is the most important class
    # add the CSF voxels, then add the GM-BG voxels, then finally the WM voxels

    WM_omean=$(mrstats -ignorezero -output mean -quiet -force ${nMNI2_inT1_ntiss_sc2T1MNI1[3]})
    
    CSF_max=$(mrstats -ignorezero -output max -quiet -force ${nMNI2_inT1_ntiss_sc2T1MNI1[0]})

    # correction for CSF in postcontrast images (seems like this is needed in general)

    if (( $(echo "${CSF_max} >= ${WM_omean}" |bc -l) )); then

        echo " Is this a postcontrast image? CSF max = ${CSF_max}, WM mean = ${WM_omean}" | tee -a ${prep_log}

    elif (( $(echo "${CSF_max} < ${WM_omean}" |bc -l) )); then

        echo " This is not a postcontrast image, CSF max = ${CSF_max}, WM mean = ${WM_omean}" | tee -a ${prep_log}

    fi

    CSF_nmean=$(mrcalc `mrstats -ignorezero -output mean -quiet -force ${nMNI2_inT1_ntiss_sc2T1MNI1[3]}` 0.01 -mult -force -quiet)

    echo " Normalized CSF tissue will be scaled to a mean of ${CSF_nmean} " | tee -a ${prep_log}

    # task_in="mrcalc -force -quiet -nthreads ${ncpu} ${MNI2_inT1_ntiss[0]} ${CSF_max} -div ${CSF_nmean} -mult ${str_pp}_nMNI2_inT1_linsc_norm_nCSF_cor.nii.gz"

    task_in="mrcalc -force -quiet -nthreads ${ncpu} ${MNI2_inT1_ntiss[0]} 0 -gt 0.2 -mult ${str_pp}_nMNI2_inT1_linsc_norm_nCSF_cor.nii.gz"

    task_exec

    nMNI2_inT1_ntiss_sc2T1MNI1[0]="${str_pp}_nMNI2_inT1_linsc_norm_nCSF_cor.nii.gz"

    GM_omean=$(mrstats -ignorezero -output mean -quiet -force ${nMNI2_inT1_ntiss_sc2T1MNI1[1]})

    WM_nmean=$(mrstats -ignorezero -output mean -quiet -force ${nMNI2_inT1_ntiss_sc2T1MNI1[3]})

    # correction for WM in low quality scans

    if (( $(echo "${GM_omean} >= ${WM_nmean}" |bc -l) )); then

        echo " Normalized WM tissue ${WM_nmean} has a lower mean than ${GM_omean} " | tee -a ${prep_log}
        echo " This could be due to low image quality " | tee -a ${prep_log}

        WM_mamean=$(mrcalc ${GM_omean} 1.75 -mult -force -quiet)

        echo " Rescaling WM tissue from ${WM_nmean} to ${WM_mamean} " | tee -a ${prep_log}

        task_in="mrcalc -force -quiet -nthreads ${ncpu} ${nMNI2_inT1_ntiss_sc2T1MNI1[3]} ${WM_nmean} -div ${WM_mamean} -mult ${str_pp}_nMNI2_inT1_linsc_norm_nWM_cor.nii.gz"

        task_exec

        nMNI2_inT1_ntiss_sc2T1MNI1[3]="${str_pp}_nMNI2_inT1_linsc_norm_nWM_cor.nii.gz"

        echo " Normalized WM tissue will be rescaled to a new mean ${WM_mamean}, which is 1.75 x times the GM mean ${GM_omean} " | tee -a ${prep_log}

    else

        # echo " Normalized WM tissue ${WM_nmean} has a lower mean than ${GM_omean} " | tee -a ${prep_log}
        # echo " Forced rescaling of WM intensity " | tee -a ${prep_log}

        # WM_mamean=$(mrcalc ${WM_nmean} 1.15 -mult -force -quiet)

        # echo " Rescaling WM tissue from a mean of ${WM_nmean} to a new mean of ${WM_mamean} " | tee -a ${prep_log}

        echo " WM tissue mean ${WM_nmean} " | tee -a ${prep_log}

        # task_in="mrcalc -force -quiet -nthreads ${ncpu} ${nMNI2_inT1_ntiss_sc2T1MNI1[3]} ${WM_nmean} -div ${WM_mamean} -mult ${str_pp}_nMNI2_inT1_linsc_norm_nWM_cor.nii.gz"

        # task_exec

        # nMNI2_inT1_ntiss_sc2T1MNI1[3]="${str_pp}_nMNI2_inT1_linsc_norm_nWM_cor.nii.gz"

        # echo " Normalized WM tissue will be rescaled to a new mean of ${WM_mamean}, which is 1.15 x times the WM old mean of ${WM_nmean} " | tee -a ${prep_log}

    fi

    task_in="ImageMath 3 ${tmp_s2T1_nCSFGMC} addtozero ${nMNI2_inT1_ntiss_sc2T1MNI1[1]} ${nMNI2_inT1_ntiss_sc2T1MNI1[0]} && ImageMath 3 \
    ${tmp_s2T1_nCSFGMCB} addtozero ${tmp_s2T1_nCSFGMC} ${nMNI2_inT1_ntiss_sc2T1MNI1[2]} && ImageMath 3 ${tmp_s2T1_nCSFGMCBWM} addtozero ${tmp_s2T1_nCSFGMCB} ${nMNI2_inT1_ntiss_sc2T1MNI1[3]}"

    task_exec

    task_in="ImageMath 3 ${tmp_s2T1_nCSFGMCBWMr} addtozero ${tmp_s2T1_nCSFGMCBWM} ${T1b_inMNI1_p_norm}"

    task_exec

    task_in="fslmaths ${tmp_s2T1_nCSFGMCBWMr} -mul ${T1b_inMNI1p_mean} ${tmp_s2T1_CSFGMCBWM}"

    task_exec

    ############

    # here we create the stitched and initial filled images
    # we also need to reconstitute the target image it seems
    # this is done in line, using fslmaths and fslstats, we divide the real image by its mean and multiply it by the mean of the synth image

    if [[ -z "${bilateral}" ]] && [[ "${t_flag}" -eq 0 ]]; then

        echo " Lesion is not bilateral ${bilateral} and no template flag is set " | tee -a ${prep_log}

        # if not bilateral, we do this using native tissue and template tissue
        
        # first we create a new healthy version of the lesioned hemisphere
        # native tissue stitched image doesnt need scaling

        task_in="fslmaths ${fT1brain_inMNI2} -mas ${L_hemi_mask} ${fT1_H_hemi}"

        task_exec

        # now we harvest the real healthy hemisphere
        
        task_in="fslmaths ${T1_brain_inMNI2} -mas ${L_hemi_mask_binv} ${T1_H_hemi}"

        task_exec

        # here we stitch the two

        task_in="fslmaths ${fT1_H_hemi} -mas ${L_hemi_mask} -add ${T1_H_hemi} ${stitched_T1_nat}"

        task_exec

        #######

        # repeat above process for template data - without filling of the donor image
        # this is not correct
        # it should be the trully scaled MNI2_in_T1
        # here we need scaling
        
        task_in="fslmaths ${tmp_s2T1_CSFGMCBWM} -mas ${L_hemi_mask} ${Temp_L_hemi}"

        task_exec

        task_in="fslmaths ${T1_brain_inMNI2} -div `mrstats -force -nthreads ${ncpu} -quiet -mask ${Lmask_binv_inMNI1_dilx2} -ignorezero -output mean ${T1_brain_inMNI2} ` \
        -mul `fslstats ${tmp_s2T1_CSFGMCBWM} -M ` -mas ${L_hemi_mask_binv} -add ${Temp_L_hemi} ${stitched_T1_temp}"

        task_exec

        #######

        stitched_T1=${stitched_T1_temp}

        # # follow same pipeline for generating a noise map

        # task_in="fslmaths ${T1_noise_inMNI1} -mas ${H_hemi_mask} ${T1_noise_H_hemi}"

        # task_exec

        # task_in="fslmaths ${fT1_noise_inMNI1} -mas ${H_hemi_mask_binv} -add ${T1_noise_H_hemi} ${stitched_noise_MNI1}"

        # task_exec

        # now we generate the initial filled map deriving the graft from the stitched_T1 using template derived tissue
        # this is to avoid midline spill over

        task_in="fslmaths ${stitched_T1} -mul ${Lmask_bin_inMNI2_s3} ${Temp_L_fill_T1}"

        task_exec


        task_in="fslmaths ${T1_brain_inMNI2} -mul ${Lmask_binv_inMNI2_s3} -add ${Temp_L_fill_T1} ${Temp_T1_filled1}"

        task_exec

        # define your initial filled map as Temp_T1_filled1
        #### I can add a move command here to simply rename this to whatever it needs to be for later on
        #### can also use a more generic name for that

        T1_filled1=${Temp_T1_filled1}

    elif [[ ! -z "${bilateral}" ]] || [[ "${t_flag}" -eq 1 ]]; then

        echo " The lesion is bilateral  bil_flag = ${bilateral} and/or template flag is set temS_flag = ${t_flag} -- " | tee -a ${prep_log}

        # similarly but no hemisphere work and no stitching
        # stitched_T1 and stitched_T1_nat are now fake ones

        stitched_T1=${tmp_s2T1_CSFGMCBWM}

        # create the initial filling graft

        task_in="fslmaths ${tmp_s2T1_CSFGMCBWM} -mul ${Lmask_bin_inMNI2_s3} ${Temp_bil_Lmask_fill1}"

        task_exec

        # Generate initial filled image

        task_in="fslmaths ${T1_brain_inMNI2} -div `mrstats -force -nthreads ${ncpu} -quiet -mask ${Lmask_binv_inMNI1_dilx2} -ignorezero -output mean ${T1_brain_inMNI2} ` \
        -mul `fslstats ${tmp_s2T1_CSFGMCBWM} -M ` -mul ${Lmask_binv_inMNI2_s3} -add ${Temp_bil_Lmask_fill1} ${Temp_T1_bilfilled1}"

        task_exec

        ######

        # if bilateral!

        T1_filled1=${Temp_T1_bilfilled1}

        stitched_T1_nat=${T1_filled1}

        echo " Since lesion is bilateral we use filled brain for Atropos1 " | tee -a ${prep_log}

    fi

    # warp the stitched to filled images, if not already done

    if [[ -z "${stitch2fill_mark}" ]] ; then

        echo "Warping stitched T1 to filled T1" | tee -a ${prep_log}

        fix_im="${T1_filled1}"

        mov_im="${stitched_T1}"

        mask=" -x ${brain_mask_inMNI1},${brain_mask_inMNI1} "

        transform="so"

        output="${T1_sti2fil_str}" 

        KUL_antsRegSyN_Def

    else

        echo " Warping stitched T1 to filled T1 already done, skipping" | tee -a ${prep_log}

    fi

    # also warp the filled T1 to the T1 in MNI1 (excluding the lesion)

    if [[ -z "${fill2MNI1mniL_mark}" ]] ; then

        echo "Warping filled T1 to T1 in MNI minL" | tee -a ${prep_log}

        fix_im="${T1_brain_inMNI1}"

        mov_im="${T1_filled1}"

        mask=" -x ${brain_mask_minL_inMNI1},${brain_mask_inMNI1} "

        transform="so"

        output="${T1fill2MNI1minL_str}" 

        KUL_antsRegSyN_Def

    else

        echo " Warping filled T1 to T1 in MNI minL already done, skipping" | tee -a ${prep_log}

    fi

    # Run Atropos 1 on the stitched ims after warping to the filled ones

    if [[ -z "${Atropos1_wf_mark}" ]] ; then

        # Run AtroposN4

        if [[ -z "${bilateral}" ]] && [[ "${t_flag}" -eq 0 ]]; then

            atropos1_brain=${T1_sti2fill_brain}

        elif [[ ! -z "${bilateral}" ]] || [[ "${t_flag}" -eq 1 ]]; then

            echo " Using the initial filled for Atropos " | tee -a ${prep_log}

            atropos1_brain=${Temp_T1_bilfilled1}

        fi

        prim_in=${atropos1_brain}

        atropos_mask="${MNI_brain_mask}"
        # this is so to avoid failures with atropos

        atropos_priors=${new_priors}

        atropos_out="${str_pp}_atropos1_"

        wt="0.3"

        mrf="[0.2,1,1,1]"

        KUL_antsAtropos

        Atropos1_str="${str_pp}_atropos1_SegmentationPosteriors?.nii.gz"

        Atropos1_posts=($(ls ${Atropos1_str}))

        echo ${Atropos1_posts[@]} | tee -a ${prep_log}

    else

        priors_str="${priors::${#priors}-9}*.nii.gz"

        # echo ${priors_str}

        priors_array=($(ls ${priors_str}))

        echo ${priors_array[@]} | tee -a ${prep_log}

        if [[ -z "${bilateral}" ]] && [[ "${t_flag}" -eq 0 ]]; then

            atropos1_brain=${T1_sti2fill_brain}

        elif [[ ! -z "${bilateral}" ]] || [[ "${t_flag}" -eq 1 ]]; then

            echo " Using the initial filled for Atropos " | tee -a ${prep_log}

            atropos1_brain=${Temp_T1_bilfilled1}

        fi

        Atropos1_str="${str_pp}_atropos1_SegmentationPosteriors?.nii.gz"

        Atropos1_posts=($(ls ${Atropos1_str}))

        echo ${Atropos1_posts[@]} | tee -a ${prep_log}

        echo " Atropos1 already done " | tee -a ${prep_log}

    fi

}

# image flipping in x with fslswapdim and fslorient 

function KUL_flip_ims {

    task_in="fslswapdim ${input} -x y z ${flip_out}"

    task_exec

    task_in="fslorient -forceradiological ${flip_out}"
    
    task_exec
    
}

# ANTs N4Atropos
# set outer iterations loop to 1 for debugging

function KUL_antsAtropos {

    task_in="antsAtroposN4.sh -d 3 -a ${prim_in} -x ${atropos_mask} -m 2 -n 6 -c 4 -y 2 -y 3 -y 4 \
    -p ${atropos_priors} -w ${wt} -r ${mrf} -o ${atropos_out} -u 1 -g 1 -k 1 -s nii.gz -z 0"

    task_exec
}


# check what kind of lesion it is

if [[ "${E_flag}" -eq 0 ]]; then

    echo
    echo "No -E flag set, treating the lesion as an intra-axial lesion, we will run VBG" >&2
    echo "No -E flag set, treating the lesion as an intra-axial lesion, we will run VBG" | tee -a ${prep_log}
    echo

    # check and report if temp flag is set

    if [[ "${t_flag}" -eq 0 ]]; then
    
        echo
        echo "Template flag not set, using native tissue for filling" >&2
        echo "Template flag not set, using native tissue for filling" | tee -a ${prep_log}
        echo
        
    elif [[ "${t_flag}" -eq 1 ]]; then

        echo
        echo " -t flag is active works best with a cooked template using KUL_VBG_cook_template.sh" >&2
        echo " -t flag is active works best with a cooked template using KUL_VBG_cook_template.sh" | tee -a ${prep_log}
        echo
        echo "Template flag is set, using native and donor tissue for filling" >&2
        echo "Template flag is set, using native and donor tissue for filling" | tee -a ${prep_log}
        echo

    fi


    # ------------------------------------------------------------------------------------------ #

    ## Start of Script

    echo " You are using VBG, please cite the following paper in your work: \
    Virtual brain grafting: Enabling whole brain parcellation in the presence of large lesions \
    Ahmed M. Radwan, Louise Emsell, Jeroen Blommaert, Andrey Zhylka, Silvia Kovacs, Tom Theys, Nico Sollmann, Patrick Dupont, Stefan Sunaert \
    medRxiv 2020.09.30.20204701; doi: https://doi.org/10.1016/j.neuroimage.2021.117731" | tee -a ${prep_log}

    echo "" | tee -a ${prep_log}

    echo " VBG started at ${start_t} " | tee -a ${prep_log}

    echo ${priors_array[@]} | tee -a ${prep_log}

    if [[ -z "${search_wf_mark1}" ]]; then

        if [[ -z "${srch_preprocp1}" ]]; then

            # inserted rescaling step 31/12/2020
            task_in="fslmaths ${prim} -thr 0.1 ${str_pp}_T1_thr.nii.gz"

            task_exec

            task_in="recon-all -i ${str_pp}_T1_thr.nii.gz -s ${subj}_temp -sd ${preproc} -openmp ${ncpu} -parallel -autorecon1"

            task_exec

            # convert to nii from the initial autorecon1 output
            task_in="mri_convert -rl ${str_pp}_T1_thr.nii.gz ${preproc}/${subj}_temp/mri/orig_nu.mgz ${str_pp}_T1_reori2std.nii.gz \
            && mri_convert -rl ${str_pp}_T1_thr.nii.gz ${preproc}/${subj}_temp/mri/brainmask.mgz ${str_pp}_brain_mask_init.nii.gz \
            && fslmaths ${str_pp}_brain_mask_init.nii.gz -bin ${str_pp}_brain_mask_FS.nii.gz"

            task_exec
            
            # && mri_convert -rl ${str_pp}_T1_thr.nii.gz ${Lmask_o} ${Lmask_FS_reori} \
            # && fslmaths ${Lmask_FS_reori} -bin ${L_mask_reori1}"            

            # initial bad BET to assist the affine reg

            # icent=$(fslstats ${str_pp}_T1_reori2std.nii.gz -C)

            # task_in="bet ${str_pp}_T1_reori2std.nii.gz ${str_pp}_T1_BET1.nii.gz -c ${icent} -f 0.45 -S -B -m -v"
            
            # task_exec

            task_in="antsRegistrationSyN.sh -d 3 -f ${MNI_T1_brain} -m ${str_pp}_brain_mask_init.nii.gz -x ${MNI_brain_mask},${str_pp}_brain_mask_FS.nii.gz \
            -o ${str_pp}_T1_reori_aff2MNI_ -t a"
    
            task_exec

            task_in="antsApplyTransforms -d 3 -i ${Lmask_o} -o ${L_mask_reori} -r ${MNI_T1_brain} -t [${str_pp}_T1_reori_aff2MNI_0GenericAffine.mat,0] -n MultiLabel \
            && antsApplyTransforms -d 3 -i ${str_pp}_T1_reori2std.nii.gz -o ${str_pp}_T1_reori_aff2MNI.nii.gz -r ${MNI_T1_brain} -t [${str_pp}_T1_reori_aff2MNI_0GenericAffine.mat,0]"

            task_exec

            task_in="DenoiseImage -d 3 -s 1 -n Gaussian -i ${str_pp}_T1_reori_aff2MNI.nii.gz -o [${str_pp}_T1_dn.nii.gz,${str_pp}_T1_noise.nii.gz] -v 1"

            task_exec

            # to avoid failures with BFC due to negative pixel values

            task_in="fslmaths ${str_pp}_T1_dn.nii.gz -thr 0.1 ${T1_pp1}"

            task_exec

        else

            echo "Reorienting, denoising, and bias correction already done, skipping " | tee -a ${prep_log}

        fi

        # Run ANTs BET, make masks

        if [[ -z "${srch_antsBET}" ]]; then

            echo " running Brain extraction " | tee -a ${prep_log}

            prim_in="${str_pp}_T1_dn_thr.nii.gz"

            output="${hdbet_str}"

            # run antsBET

            KUL_antsBETp

            task_in="fslmaths ${clean_mask_nat} -s 2 -thr 0.5 -save ${BET_mask_s2} -binv -fillh -s 2 -thr 0.2 \
            -sub ${BET_mask_s2} -thr 0 ${BET_mask_binvs2}"

            task_exec

            task_in="fslmaths ${T1_pp1} -mul ${BET_mask_binvs2} ${T1_skull}"

            task_exec

        else

            echo " Brain extraction already done, skipping " | tee -a ${prep_log}

            echo "${T1_brain_clean}" | tee -a ${prep_log}

            echo "${clean_mask_nat}" | tee -a ${prep_log}

            echo " ANTsBET already run, skipping " | tee -a ${prep_log}

        fi

        # exit 2

        # exit if theres a problem with brain extraction

        srch_antsBET2=($(find ${preproc} -type f | grep "${T1_brain_clean}"));

        if [[ -z "${srch_antsBET2}" ]]; then

            echo " Brain extraction not successful, please see logs exiting"
            echo " Brain extraction not successful, please see logs exiting" | tee -a ${prep_log}
            exit 2

        else

            echo " Brain extraction successful, carry on" | tee -a ${prep_log}

        fi

        # run KUL_lesion_magic1
        # this creates a bin, binv, & bm_minL 

        if [[ -z "${sch_brnmsk_minL}" ]]; then

            if [[ "${E_flag}" -eq 0 ]]; then

                # echo " Lesion mask is already in T1 space " | tee -a ${prep_log}

                echo " Intra-axial lesion running VBG Lmask_pt1 workflow and subsequent steps" | tee -a ${prep_log}

                # start by smoothing and thring the mask

                task_in="fslmaths ${L_mask_reori} -s 2 -thr 0.2 -bin -save ${Lmask_bin} -binv ${Lmask_in_T1_binv}"

                task_exec

                echo " Copying Lmask_bin_s2 file to Lmask_in_T1_bin " | tee -a ${prep_log}

                cp ${Lmask_bin} ${Lmask_in_T1_bin}

                # subtract lesion from brain mask

                task_in="fslmaths ${clean_mask_nat} -mas ${Lmask_in_T1_binv} -mas ${clean_mask_nat} ${brain_mask_minL}"

                task_exec

            else

                echo " Extra-axial lesion running simplified VBG Lmask_pt1 workflow, FS and subsequent steps" | tee -a ${prep_log}

                task_in="fslmaths ${L_mask_reori} -binv ${L_O_binv}"

                task_exec

            
            fi

        else

            echo "${brain_mask_minL} already created " | tee -a ${prep_log}

        fi

        # carry on
        # here we do the first T1 warp to template (default antsRegSyN 3 stage)

        if [[ -z "${T1brain2MNI1}" ]] ; then
        
            mov_im="${T1_brain_clean}"

            fix_im="${MNI_T1_brain}"

            mask=" -x ${MNI_brain_mask},${brain_mask_minL} "

            transform="s"

            output="${T1_brMNI1_str}"

            KUL_antsRegSyN_Def

        else

            echo "${T1_brain_inMNI1} already created " | tee -a ${prep_log}

        fi

        # Apply warps to brain_mask, noise, L_mask and make BM_minL

        task_in="antsApplyTransforms -d 3 -i ${clean_mask_nat} -o ${brain_mask_inMNI1} -r ${MNI_T1_brain} -t ${T1_brMNI1_str}1Warp.nii.gz -t [${T1_brMNI1_str}0GenericAffine.mat,0] -n MultiLabel"

        task_exec

        # task_in="WarpImageMultiTransform 3 ${str_pp}_T1_noise.nii.gz ${T1_noise_inMNI1} -R ${MNI_T1_brain} ${T1_brMNI1_str}1Warp.nii.gz ${T1_brMNI1_str}0GenericAffine.mat --use-NN"

        # task_exec

        task_in="antsApplyTransforms -d 3 -i ${Lmask_in_T1_bin} -o ${str_pp}_Lmask_rsMNI1.nii.gz -r ${MNI_T1_brain} -t ${T1_brMNI1_str}1Warp.nii.gz -t [${T1_brMNI1_str}0GenericAffine.mat,0] -n MultiLabel"

        task_exec
        
        task_in="fslmaths ${str_pp}_Lmask_rsMNI1.nii.gz -bin -mas ${brain_mask_inMNI1} -save ${Lmask_bin_inMNI1} -binv -mas ${brain_mask_inMNI1} ${brain_mask_minL_inMNI1}"

        task_exec


    else


        echo " First part already done, skipping. " | tee -a ${prep_log}

    fi

    # Flip the images in MNI1

    if [[ -z "${search_wf_mark2}" ]]; then

        # do it for the T1s

        input="${T1_brain_inMNI1}"

        flip_out="${fT1brain_inMNI1}"

        KUL_flip_ims

        unset input flip_out

        input="${brain_mask_minL_inMNI1}"

        flip_out="${fbrain_mask_minL_inMNI1}"

        KUL_flip_ims

        unset input flip_out

    else

        echo " flipped images already created, skipping " | tee -a ${prep_log}

    fi

    # Second deformation, warp to template a second time (1 stage SyN)

    if [[ -z "${T1brain2MNI2}" ]]; then

        fix_im="${MNI_T1_brain}"

        mov_im="${T1_brain_inMNI1}"

        transform="s"

        output="${T1_brMNI2_str}"

        mask=" -x ${MNI_brain_mask},${brain_mask_minL_inMNI1} "

        KUL_antsRegSyN_Def

    else

        echo " Second T1_brain 2 MNI already done, skipping " | tee -a ${prep_log}

    fi

    # Warp flipped brain to template (3 stage)

    if [[ -z "${fT1_brain_2MNI2}" ]]; then

        fix_im="${MNI_T1_brain}"

        mov_im="${fT1brain_inMNI1}"

        transform="s"

        output="${fT1_brMNI2_str}"

        mask=" -x ${MNI_brain_mask},${fbrain_mask_minL_inMNI1} "

        echo " now making fT1brain2MNI2 " | tee -a ${prep_log}

        KUL_antsRegSyN_Def

    else

        echo " Second fT1_brain 2 MNI already done, skipping " | tee -a ${prep_log}

    fi

    # Atropos2 runs on the T1 brain in MNI1 space (after first deformation)
    # after this runs, cleanup and apply inverse warps to native space

    if [[ -z "${Atropos2_wf_mark}" ]]; then

        unset atropos_out prim_in atropos_mask atropos_priors wt mrf

        echo " using MNI priors for segmentation "  | tee -a ${prep_log}

        echo ${priors_array[@]} | tee -a ${prep_log}

        prim_in=${T1_brain_inMNI1}

        atropos_mask="${brain_mask_minL_inMNI1}"

        echo " using brain mask minL in MNI1 "  | tee -a ${prep_log}
        # this is so to avoid failures with atropos

        atropos_priors=${new_priors}

        atropos_out="${str_pp}_atropos2_"

        wt="0.1"

        mrf="[0.1,1,1,1]"

        KUL_antsAtropos

        Atropos2_str="${str_pp}_atropos2_SegmentationPosteriors?.nii.gz"

        Atropos2_posts=($(ls ${Atropos2_str}))

        echo ${Atropos2_posts[@]} | tee -a ${prep_log}

        unset atropos_out prim_in atropos_mask atropos_priors wt mrf

    else

        Atropos2_str="${str_pp}_atropos2_SegmentationPosteriors?.nii.gz"

        Atropos2_posts=($(ls ${Atropos2_str}))

        echo ${Atropos2_posts[@]} | tee -a ${prep_log}

        echo " Atropos2 segmentation already finished, skipping. " | tee -a ${prep_log}

    fi

    # second part handling lesion masks and Atropos run
    # this function has builtin processing control points

    if [[ ! -f "${str_pp}_atropos1_Segmentation.nii.gz" ]]; then

        echo " Starting KUL lesion magic part 2 " | tee -a ${prep_log}

        KUL_Lmask_part2

        echo " Finished KUL lesion magic part 2 " | tee -a ${prep_log}

    else

        task_in="fslmaths ${MNI_lwr} -mas ${Lmask_bin_inMNI1} ${lesion_left_overlap}"

        task_exec

        task_in="fslmaths ${MNI_rwr} -mas ${Lmask_bin_inMNI1} ${lesion_right_overlap}"

        task_exec

        Lmask_tot_v=$(mrstats -force -nthreads ${ncpu} ${Lmask_bin_inMNI1} -output count -quiet -ignorezero)

        overlap_left=$(mrstats -force -nthreads ${ncpu} ${lesion_left_overlap} -output count -quiet -ignorezero)

        overlap_right=$(mrstats -force -nthreads ${ncpu} ${lesion_right_overlap} -output count -quiet -ignorezero)

        T1b_inMNI1p_mean=$(fslstats ${T1b_inMNI1_punched} -M)

        L_ovLt_2_total=$(echo ${overlap_left}*100/${Lmask_tot_v} | bc)

        L_ovRt_2_total=$(echo ${overlap_right}*100/${Lmask_tot_v} | bc)

        echo " total lesion vox count ${Lmask_tot_v}" | tee -a ${prep_log}

        echo " ov_left is ${overlap_left}" | tee -a ${prep_log}

        echo " ov right is ${overlap_right}" | tee -a ${prep_log}

        echo " ov_Lt to total is ${L_ovLt_2_total}" | tee -a ${prep_log}

        echo " ov_Rt to total is ${L_ovRt_2_total}" | tee -a ${prep_log}

        # we set a hard-coded threshold of 65, if unilat. then native heatlhy hemi is used
        # if bilateral by more than 35, template brain is used
        # # this needs to be modified, also need to include simple lesion per hemisphere overlap with percent to total hemi volume
        # this will enable us to use template or simple filling and derive mean values per tissue class form another source (as we are currently using the original images).
        # AR 09/02/2020
        # here we also need to make unilateral L masks, masked by hemi mask to overcome midline issue
        
        if [[ "${L_ovLt_2_total}" -gt 65 ]]; then

            echo ${L_ovLt_2_total} | tee -a ${prep_log}

            echo " This patient has a left sided or predominantly left sided lesion " | tee -a ${prep_log}

            echo "${L_hemi_mask}" | tee -a ${prep_log}

            echo "${H_hemi_mask}" | tee -a ${prep_log}

            stitched_T1=${tmp_s2T1_CSFGMCBWM}

            stitched_T1=${stitched_T1_temp}

        elif [[ "${L_ovRt_2_total}" -gt 65 ]]; then

            echo ${L_ovRt_2_total} | tee -a ${prep_log}

            echo " This patient has a right sided or predominantly right sided lesion " | tee -a ${prep_log}

            echo "${L_hemi_mask}" | tee -a ${prep_log}

            echo "${H_hemi_mask}" | tee -a ${prep_log}

            T1_filled1=${Temp_T1_filled1}
            
            stitched_T1=${tmp_s2T1_CSFGMCBWM}
            
        else 

            bilateral=1

            stitched_T1=${stitched_T1_temp}

            T1_filled1=${Temp_T1_bilfilled1}

            stitched_T1_nat=${T1_filled1}

            echo " This is a bilateral lesion with ${L_ovLt_2_total} left side and ${L_ovRt_2_total} right side, using Template T1 to derive lesion fill patch. "  | tee -a ${prep_log}

            echo " note Atropos1 will use the filled images instead of the stitched ones "  | tee -a ${prep_log}

        fi

        if [[ -z "${bilateral}" ]] && [[ "${t_flag}" -eq 0 ]]; then

            atropos1_brain=${T1_sti2fill_brain}

        elif [[ ! -z "${bilateral}" ]] || [[ "${t_flag}" -eq 1 ]]; then

            echo " Using the initial filled for Atropos " | tee -a ${prep_log}

            atropos1_brain=${Temp_T1_bilfilled1}

        fi

        Atropos1_str="${str_pp}_atropos1_SegmentationPosteriors?.nii.gz"

        Atropos1_posts=($(ls ${Atropos1_str}))

        echo ${Atropos1_posts[@]} | tee -a ${prep_log}

        echo " Lesion magic part 2 already finished, skipping " | tee -a ${prep_log}


    fi

    ##### 

    # Now we warp back to MNI brain in native space
    # which will be needed after Atropos1
    # this can be replaced by a different kind of reg no ? or we simply apply the inverse warps!
    # just to show the initial filled result (initially for diagnostic purposes)

    if [[ -z "${srch_bk2anat1_mark}" ]]; then

        # fix_im="${T1_brMNI1_str}InverseWarped.nii.gz"

        fix_im="${T1_filled1}"

        mov_im="${T1_brain_clean}"

        mask=" -x ${brain_mask_inMNI1},${brain_mask_minL} "

        transform="s"

        output="${T1_bk2nat1_str}"

        KUL_antsRegSyN_Def

        task_in="cp ${T1_bk2nat1_str}InverseWarped.nii.gz ${str_op}_T1_initial_filled_brain.nii.gz"

        task_exec    

        # the bk2anat1 step is for the outputs of Atropos1 mainly.

    else

        echo "First step Warping images back to anat already done, skipping " | tee -a ${prep_log}

    fi

    # we will need fslmaths to make lesion_fill2 from the segmentations
    # fslmaths Atropos2_posteriors -add lesion_fill2
    # fslstats -m
    # fslmaths to binarize each tpm then -mul the mean intensity of that tissue type
    # should mask out the voxels of each tpm from the resulting image before inserting it
    # finally fslmaths -add noise and -mul bias

    echo " Starting KUL lesion magic part 3 " | tee -a ${prep_log}

    if [[ -z "${srch_make_images}" ]]; then 

        # warping nat filled (in case of unilat. lesion and place holder for initial filled)
        # will be used to fill holes in synth image

        # task_in="WarpImageMultiTransform 3 ${stitched_T1_nat} ${stitched_T1_nat_innat} -R ${T1_brain_clean} \
        # -i ${T1_bk2nat1_str}0GenericAffine.mat ${T1_bk2nat1_str}1InverseWarp.nii.gz && WarpImageMultiTransform 3 ${stitched_T1_temp} \
        # ${stitched_T1_temp_innat} -R ${T1_brain_clean} -i ${T1_bk2nat1_str}0GenericAffine.mat ${T1_bk2nat1_str}1InverseWarp.nii.gz"

        task_in="antsApplyTransforms -d 3 -i ${stitched_T1} -o ${stitched_T1_innat} -r ${T1_brain_clean} -t [${T1_bk2nat1_str}0GenericAffine.mat,1] -t ${T1_bk2nat1_str}1InverseWarp.nii.gz"

        task_exec

        # Create the segmentation image lesion fill

        task_in="fslmaths ${str_pp}_atropos1_Segmentation.nii.gz -mas ${Lmask_bin_inMNI1_dilx2} ${Lfill_segm_im}"

        task_exec

        # Make a hole in the real segmentation image and fill it

        task_in="fslmaths ${str_pp}_atropos2_Segmentation.nii.gz -mas ${Lmask_binv_inMNI1_dilx2} -add ${Lfill_segm_im} ${atropos2_segm_im_filled}"

        task_exec

        # Bring the Atropos2_segm_im to native space
        # for diagnostic purposes

        task_in="antsApplyTransforms -d 3 -i ${atropos2_segm_im_filled} -o ${atropos2_segm_im_filled_nat} -r ${T1_brain_clean} -t [${T1_bk2nat1_str}0GenericAffine.mat,1] -t ${T1_bk2nat1_str}1InverseWarp.nii.gz -n MultiLabel"

        task_exec

        # Making the cleaned segmentation images here

        echo ${tissues[@]} | tee -a ${prep_log}

        # T1b_inMNI1_p_norm is generated before 

        # Match the intensities of the normalized tissue components from each image

        image_out="${str_pp}_donor_T1_native.nii.gz"

        task_in="antsApplyTransforms -d 3 -i ${T1_sti2fill_brain} -o ${image_out} -r ${T1_brain_clean} -t [${T1_brMNI1_str}0GenericAffine.mat,1] -t ${T1_brMNI1_str}1InverseWarp.nii.gz"

        task_exec

        unset image_in image_out
        
        task_in="ImageMath 3 ${str_pp}_donor_T1_native_S.nii.gz Sharpen ${str_pp}_donor_T1_native.nii.gz"

        task_exec

        task_in="fslmaths ${str_pp}_donor_T1_native_S.nii.gz -div `mrstats -force -nthreads ${ncpu} -quiet -mask ${clean_mask_nat} -ignorezero -output mean ${str_pp}_donor_T1_native_S.nii.gz ` \
        -mul `mrstats -mask ${brain_mask_minL} -force -nthreads ${ncpu} -quiet -ignorezero -output mean ${T1_brain_clean} ` \
        -save ${str_op}_donor_brain.nii.gz -mul ${Lmask_bin_s3} ${T1_fin_Lfill_1}"

        task_exec
        
        # make the final outputs
    
        task_in="fslmaths ${T1_brain_clean} -mul ${Lmask_binv_s3} -add ${T1_fin_Lfill_1} -thr 0 -save ${T1_nat_filled_out_1} -mul ${BET_mask_s2} \
        -add ${T1_skull} -thr 0 ${T1_nat_fout_wskull_1} && ImageMath 3 ${T1_nat_filled_out_2} HistogramMatch ${T1_nat_filled_out_1} ${T1_pp1} \
        && ImageMath 3 ${T1_nat_fout_wskull_2} HistogramMatch ${T1_nat_fout_wskull_1} ${T1_pp1} \
        && fslmaths ${T1_nat_fout_wskull_2} -mul ${clean_mask_nat} ${T1_nat_filled_out_2}"

        task_exec

        task_in="ImageMath 3 ${Lmask_bin_s3_flat} FlattenImage ${Lmask_bin_s3} 2 && fslmaths ${T1_nat_fout_wskull_2} -mul ${Lmask_bin_s3_flat} ${T1_fin_Lfill_2} \
        && fslmaths ${clean_mask_nat} -mul 0 -add 1 -sub ${Lmask_bin_s3} ${Lmask_binv_s3_nobrain}"

        task_exec
        
        task_in="antsApplyTransforms -d 3 -i ${Lmask_binv_s3_nobrain} -o ${Lmask_binv_s3_n_ori} -r ${str_pp}_brain_mask_init.nii.gz -t [${str_pp}_T1_reori_aff2MNI_0GenericAffine.mat,1] \
        && antsApplyTransforms -d 3 -i ${T1_fin_Lfill_2} -o ${T1_fin_Lfill_n_ori} -r ${str_pp}_brain_mask_init.nii.gz -t [${str_pp}_T1_reori_aff2MNI_0GenericAffine.mat,1] \
        && antsApplyTransforms -d 3 -i ${clean_mask_nat} -o ${T1_BM_4_FS} -r ${str_pp}_brain_mask_init.nii.gz -t [${str_pp}_T1_reori_aff2MNI_0GenericAffine.mat,1] -n MultiLabel \
        && fslmaths ${str_pp}_T1_reori2std.nii.gz -mul ${Lmask_binv_s3_n_ori} -add ${T1_fin_Lfill_n_ori} -thr 0 -save ${T1_4_FS} -mul ${T1_BM_4_FS} \
        ${T1_Brain_4_FS} && mri_convert -i ${T1_4_FS} -o ${T1_4_parc} --conform"

        task_exec
        
        # task_in="convert_xfm -omat ${T1_reori_mat_inv} -inverse ${T1_reori_mat} && sleep 5 \
        # && flirt -in ${Lmask_binv_s3_nobrain} -out ${Lmask_binv_s3_n_ori} -ref ${T1_orig} -applyxfm -init ${T1_reori_mat_inv} \
        # && sleep 5 && flirt -in ${T1_fin_Lfill_2} -out ${T1_fin_Lfill_n_ori} -ref ${T1_orig} -applyxfm -init ${T1_reori_mat_inv} \
        # && sleep 5 && flirt -in ${clean_mask_nat} -out ${T1_BM_4_FS} -ref ${T1_orig} -applyxfm -init ${T1_reori_mat_inv} \
        # && fslmaths ${str_pp}_T1_thr.nii.gz -mul ${Lmask_binv_s3_n_ori} -add ${T1_fin_Lfill_n_ori} -thr 0 -save ${T1_4_FS} -mul ${T1_BM_4_FS} \
        # ${T1_Brain_4_FS}"

        # task_exec

    else

        echo " Making fake healthy images done, skipping. " | tee -a ${prep_log}

    fi


    unset i
    
else

    echo
    echo "You have set the -E flag, indicating an extra-axial lesion" 
    echo "You have set the -E flag, indicating an extra-axial lesion" | tee -a ${prep_log}
    echo "The lesion patch is filled with 0s only, recon-all should be able to run, if it fails try without -E" 
    echo "The lesion patch is filled with 0s only, recon-all should be able to run, if it fails try without -E" | tee -a ${prep_log}
    echo

    if [[ -z "${srch_preprocp1}" ]]; then

            input="${prim}"

            output1="${str_pp}_T1_reori2std.nii.gz"

            # task_in="fslreorient2std -m ${T1_reori_mat} ${input} ${output1}"

            task_in="fslreorient2std ${input} ${output1} && fslreorient2std ${input} >> ${T1_reori_mat}"

            task_exec

        if [[ -z "${srch_antsBET}" ]]; then

            echo " running Brain extraction " | tee -a ${prep_log}

            prim_in="${str_pp}_T1_reori2std.nii.gz"

            output="${hdbet_str}"

            # run antsBET

            KUL_antsBETp


        else

            echo " Brain extraction already done, skipping " | tee -a ${prep_log}

            echo "${T1_brain_clean}" | tee -a ${prep_log}

            echo "${clean_mask_nat}" | tee -a ${prep_log}

            echo " ANTsBET already run, skipping " | tee -a ${prep_log}

        fi

        # run KUL_lesion_magic1
        # this creates a bin, binv, & bm_minL 

        if [[ -z "${sch_brnmsk_minL}" ]]; then

            KUL_Lmask_part1

        else

            echo "${brain_mask_minL} already created " | tee -a ${prep_log}

        fi

    else

        echo "Reorienting, brain extraction, and VBG pt1 already done, skipping " | tee -a ${prep_log}
    
        
    fi

    # here we fill the lesion mask with 0 and save it where FS recon-all expects it to be

    task_in="fslmaths ${str_pp}_T1_reori2std.nii.gz -mul ${L_O_binv} ${T1_nat_fout_wskull_2}"

    task_exec

    task_in="antsApplyTransforms -d 3 -i ${Lmask_binv_s3_nobrain} -o ${Lmask_binv_s3_n_ori} -r ${str_pp}_brain_mask_init.nii.gz -t [${str_pp}_T1_reori_aff2MNI_0GenericAffine.mat,1] \
    && antsApplyTransforms -d 3 -i ${T1_fin_Lfill_2} -o ${T1_fin_Lfill_n_ori} -r ${str_pp}_brain_mask_init.nii.gz -t [${str_pp}_T1_reori_aff2MNI_0GenericAffine.mat,1] \
    && antsApplyTransforms -d 3 -i ${clean_mask_nat} -o ${T1_BM_4_FS} -r ${str_pp}_brain_mask_init.nii.gz -t [${str_pp}_T1_reori_aff2MNI_0GenericAffine.mat,1] -n MultiLabel \
    && fslmaths ${str_pp}_T1_reori2std.nii.gz -mul ${Lmask_binv_s3_n_ori} -add ${T1_fin_Lfill_n_ori} -thr 0 -save ${T1_4_FS} -mul ${T1_BM_4_FS} \
    ${T1_Brain_4_FS} && mri_convert -i ${T1_4_FS} -o ${T1_4_parc} --conform"

    task_exec

fi



# now we need to try all the above steps, and debug, then program a function for the lesion patch filling
# then add in the recon-all step
# and add in again the overlap calculator and report parts.

# classic_FS

# for recon-all

if [[ "${P_flag}" -eq 1 ]] ; then

    if [[ "${parc_F}" -eq 1 ]] ; then

        echo
        echo "Fresurfer flag is set, now starting FS recon-all based part of VBG" >&2
        echo "Fresurfer flag is set, now starting FS recon-all based part of VBG" | tee -a ${prep_log}
        echo

        if [[ "$bids_flag" -eq 1 ]] && [[ "$o_flag" -eq 0 ]]; then

            fs_output="${cwd}/BIDS/derivatives/freesurfer/sub-${subj}"

        else

            fs_output="${str_op}_FS_output/sub-${subj}"

        fi

        recall_scripts="${fs_output}/${subj}/scripts"

        search_wf_mark4=($(find ${recall_scripts} -type f 2> /dev/null | grep recon-all.done));

        FS_brain="${fs_output}/${subj}/mri/brainmask.mgz"

        new_brain="${str_pp}_T1_Brain_4FS.mgz"

        if [[ -z "${search_wf_mark4}" ]]; then

            task_in="mkdir -p ${fs_output} >/dev/null 2>&1"

            task_exec

            # Run recon-all and convert the real T1 to .mgz for display
            # running with -noskulltrip and using brain only inputs
            # for recon-all
            # if we can run up to skull strip, break, fix with hd-bet result then continue it would be much better
            # if we can switch to fast-surf, would be great also
            # another possiblity is using recon-all -skullstrip -clean-bm -gcut -subjid <subject name>

            task_in="recon-all -i ${T1_4_parc} -s ${subj} -sd ${fs_output} -openmp ${ncpu} -parallel -autorecon1"

            task_exec

            task_in="mri_convert -rl ${fs_output}/${subj}/mri/brainmask.mgz ${T1_BM_4_FS} ${clean_BM_mgz}"

            task_exec

            task_in="mri_mask ${FS_brain} ${T1_BM_4_FS} ${new_brain} && mv ${new_brain} ${fs_output}/${subj}/mri/brainmask.mgz && cp \
            ${fs_output}/${subj}/mri/brainmask.mgz ${fs_output}/${subj}/mri/brainmask.auto.mgz"

            task_exec

            task_in="recon-all -s ${subj} -sd ${fs_output} -openmp ${ncpu} -parallel -all -noskullstrip"

            task_exec

            task_in="mri_convert -rl ${fs_output}/${subj}/mri/brain.mgz ${T1_brain_clean} ${fs_output}/${subj}/mri/real_T1.mgz"

            task_exec

            task_in="mri_convert -rl ${fs_output}/${subj}/mri/brain.mgz -rt nearest ${Lmask_o} ${fs_output}/${subj}/mri/Lmask_T1_bin.mgz"

            task_exec

            fs_parc_mgz="${fs_output}/${subj}/mri/aparc+aseg.mgz"

        else

            echo " recon-all already done, skipping. "
            echo " recon-all already done, skipping. "  | tee -a ${prep_log}
            
            fs_parc_mgz="${fs_output}/${subj}/mri/aparc+aseg.mgz"
            
        fi

    elif [[ "${parc_F}" -eq 2 ]] ; then

        echo
        echo "FastSurfer flag is set, now starting FaSu recon-all based part of VBG" >&2
        echo "FastSurfer flag is set, now starting FaSu recon-all based part of VBG" | tee -a ${prep_log}
        echo

        if [[ "$bids_flag" -eq 1 ]] && [[ "$o_flag" -eq 0 ]]; then

            fs_output="${cwd}/BIDS/derivatives/fastsurfer/sub-${subj}"

        else

            fs_output="${str_op}_FaSu_output/sub-${subj}"

        fi

        recall_scripts="${fs_output}/${subj}/scripts"

        search_wf_mark4=($(find ${recall_scripts} -type f 2> /dev/null | grep recon-all.done));

        FS_brain="${fs_output}/${subj}/mri/brainmask.mgz"

        new_brain="${str_pp}_T1_Brain_4FS.mgz"

        if [[ ! -f "${fs_output}/${subj}/label/rh.aparc.annot" ]]; then
        
            if [[ -z "${search_wf_mark4}" ]]; then

                task_in="mkdir -p ${fs_output} >/dev/null 2>&1"

                task_exec

                # Now we need to figure out how to run FaSu
                # this is really a beta option
                # so set it according to your arch ?
                # a general solution is the fastsurfer cpu docker file

                # search for FaSu native install first
                FaSu_loc=$(which run_fastsurfer.sh)
                nvd_cu=$(nvcc --version)
                user_id_str=$(id -u $(whoami))

                if [[ ! -z ${FaSu_loc} ]]; then

                    if [[ -z ${nvd_cu} ]]; then

                        FaSu_cpu=" --no_cuda "

                    else

                        FaSu_cpu=""

                    fi

                    # it's a good idea to run autorecon1 first anyway
                    # then use the orig from that to feed to FaSu

                    task_in="run_fastsurfer.sh --t1 ${T1_4_parc} \
                    --sid ${subj} --sd ${fs_output} --parallel --fsaparc --threads ${ncpu} \
                    --fs_license $FREESURFER_HOME/license.txt --py python ${FaSu_cpu}"

                    task_exec

                    task_in="mri_convert -rl ${fs_output}/${subj}/mri/brain.mgz ${T1_brain_clean} ${fs_output}/${subj}/mri/real_T1.mgz"

                    task_exec

                    task_in="mri_convert -rl ${fs_output}/${subj}/mri/brain.mgz -rt nearest ${Lmask_o} ${fs_output}/${subj}/mri/Lmask_T1_bin.mgz"

                    task_exec

                else

                    # it's a good idea to run autorecon1 first anyway
                    # then use the orig from that to feed to FaSu

                    echo "Local FastSurfer not found, switching to Docker version" | tee -a ${prep_log}
                    T1_4_FaSu=$(basename ${T1_4_parc})

                    if [[ ! -z ${nvd_cu} ]]; then

                        FaSu_v="gpu"

                    else

                        FaSu_v="cpu"

                    fi

                    task_in="docker run -v ${output_d}:/data -v ${fs_output}:/output \
                    -v $FREESURFER_HOME:/fs60 --rm --user ${user_id_str} fastsurfer:${FaSu_v} \
                    --fs_license /fs60/license.txt --sid ${subj} \
                    --sd /data/sub-${subj}${ses_long}_FS_output/sub-${subj} --t1 /data/${T1_4_FaSu} \
                    --parallel --fsaparc --threads ${ncpu}"

                    task_exec

                    task_in="mri_convert -rl ${fs_output}/${subj}/mri/brain.mgz ${T1_brain_clean} ${fs_output}/${subj}/mri/real_T1.mgz"

                    task_exec

                    task_in="mri_convert -rl ${fs_output}/${subj}/mri/brain.mgz -rt nearest ${Lmask_o} ${fs_output}/${subj}/mri/Lmask_T1_bin.mgz"

                    task_exec


                fi
                # check if cuda is available
                # if so do FaSu natively with CUDA
                # if not found then search for docker
                # check if cuda is available
                # if both found then use 
                # docker with cuda
                # if cuda not found then docker cpu

            fi

            # need to find fs wd and go to the level of license

            fs_parc_mgz="${fs_output}/${subj}/mri/aparc.mapped+aseg.mgz"

            cp ${fs_output}/${subj}/label/lh.aparc.mapped.annot ${fs_output}/${subj}/label/lh.aparc.annot

            cp ${fs_output}/${subj}/label/rh.aparc.mapped.annot ${fs_output}/${subj}/label/rh.aparc.annot
        
        else

            echo " recon-all already done, skipping. "
            echo " recon-all already done, skipping. "  | tee -a ${prep_log}
            fs_parc_mgz="${fs_output}/${subj}/mri/aparc.mapped+aseg.mgz"

            #cp ${fs_output}/${subj}/label/lh.aparc.mapped.annot ${fs_output}/${subj}/label/lh.aparc.annot

            #cp ${fs_output}/${subj}/label/rh.aparc.mapped.annot ${fs_output}/${subj}/label/rh.aparc.annot
            
        fi


    fi

    # ## After recon-all is finished we need to calculate percent lesion/lobe overlap
    # # need to make labels array

    lesion_lobes_report="${fs_output}/percent_lobes_lesion_overlap_report.txt"

    task_in="touch ${lesion_lobes_report}"

    task_exec

    echo " Percent overlap between lesion and each lobe " | tee -a $lesion_lobes_report

    echo " each lobe mask voxel count and volume in cmm is reported " | tee -a $lesion_lobes_report

    echo " overlap in voxels and volume cmm are reported " | tee -a $lesion_lobes_report

    # these labels, wm and gm values are used later for the reporting

    # double checking: RT_Frontal, LT_Frontal, RT_Temporal, LT_Temporal 

    declare -a labels=("RT_Frontal"  "LT_Frontal"  "RT_Temporal"  "LT_Temporal"  "RT_Parietal"  "LT_Parietal" \
    "RT_Occipital"  "LT_Occipital"  "RT_Cingulate"  "LT_Cingulate"  "RT_Insula"  "LT_Insula"  "RT_Putamen"  "LT_Putamen" \
    "RT_Caudate"  "LT_Caudate"  "RT_Thalamus"  "LT_Thalamus" "RT_Pallidum"  "LT_Pallidum"  "RT_Accumbens"  "LT_Accumbens"  "RT_Amygdala"  "LT_Amygdala" \
    "RT_Hippocampus"  "LT_Hippocampus"  "RT_PWM"  "LT_PWM");

    declare -a wm=("4001"  "3001"  "4005"  "3005"  "4006"  "3006" \
    "4004"  "3004"  "4003"  "3003"  "4007"  "3007" "0"  "0" \
    "0"  "0"  "0"  "0"  "0"  "0"  "0"  "0"  "0"  "0" \
    "0"  "0"  "5002"  "5001");

    declare -a gm=("2001"  "1001"  "2005"  "1005"  "2006"  "1006" \
    "2004"  "1004"  "2003"  "1003"  "2007"  "1007" "51"  "12" \
    "50"  "11"  "49"  "10"  "52"  "13"  "58"  "26" "54"  "18" \
    "53"  "17"  "0"  "0");

    fs_lobes_mgz="${fs_output}/${subj}/mri/lobes_ctx_wm_fs.mgz"

    fs_parc_nii="${str_op}_aparc+aseg.nii.gz"

    fs_parc_minL_nii="${str_op}_aparc+aseg_minL.nii.gz"

    fs_lobes_nii="${str_op}_lobes_ctx_wm_fs.nii.gz"

    fs_lobes_minL_nii="${str_op}_lobes_ctx_wm_fs_minL.nii.gz"

    fs_parc_plusL_nii="${str_op}_aparc+aseg+Lesion.nii.gz"

    fs_lobes_plusL_nii="${str_op}_lobes_ctx_wm_fs+Lesion.nii.gz"

    fs_parc_nii_LC="${str_op}_aparc_LC.nii.gz"

    fs_parc_plusL_nii_LC="${str_op}_aparc+Lesion_LC.nii.gz"

    fs_parc_minL_nii_LC="${str_op}_aparc_minL_LC.nii.gz"

    labelslength=${#labels[@]}

    wmslength=${#wm[@]}

    gmslength=${#gm[@]}

    fs_lobes_mark=${fs_lobes_nii}

    search_wf_mark5=($(find ${output_d} -type f | grep lobes_ctx_wm_fs+Lesion.nii));

    if [[ -z "$search_wf_mark5" ]]; then

        # quick sanity check

        if [[ "${labelslength}" -eq "${wmslength}" ]] && [[ "${gmslength}" -eq "${wmslength}" ]]; then

            echo "we are doing okay captain! ${labelslength} ${wmslength} ${gmslength}" | tee -a ${prep_log}

        else

            echo "we have a problem captain! ${labelslength} ${wmslength} ${gmslength}" | tee -a ${prep_log}
            
            exit 2

        fi

        # this approach apparently screws up the labels order, so i need to use annotation2label and mergelabels instead.

        task_in="mri_annotation2label --subject ${subj} --sd ${fs_output} --hemi rh --lobesStrict ${fs_output}/${subj}/label/rh.lobesStrict"

        task_exec

        task_in="mri_annotation2label --subject ${subj} --sd ${fs_output} --hemi lh --lobesStrict ${fs_output}/${subj}/label/lh.lobesStrict"

        task_exec

        task_in="mri_aparc2aseg --s ${subj} --sd ${fs_output} --labelwm --hypo-as-wm --rip-unknown --volmask --annot lobesStrict --o ${fs_lobes_mgz}"

        task_exec

        task_in="mri_convert -rl ${T1_4_FS} -rt nearest ${fs_lobes_mgz} ${fs_lobes_nii}"
        
        task_exec

        task_in="mri_convert -rl ${T1_4_FS} -rt nearest ${fs_parc_mgz} ${fs_parc_nii}"
        
        task_exec

        # here we want to add a loop looking at lesion mask volume
        l_vol=($(fslstats ${Lmask_o} -V))

        # echo "this lesion is not larger than 10 ml, we will not erode it"

        task_in="fslmaths ${Lmask_o} -binv -mul ${T1_BM_4_FS} -bin ${bmc_minL_true}"

        task_exec

        task_in="maskfilter -force -nthreads ${ncpu} ${bmc_minL_true} connect - -connectivity -largest | mrcalc - 0.1 -gt ${bmc_minL_conn} -force -nthreads ${ncpu} -quiet"

        task_exec

        task_in="fslmaths ${Lmask_o} -bin -mul 99 ${L_mask_reori_scaled}"

        task_exec

        task_in="fslmaths ${fs_parc_nii} -mas ${bmc_minL_conn} ${fs_parc_minL_nii} && ImageMath 3 ${fs_parc_plusL_nii} \
        addtozero ${fs_parc_minL_nii} ${L_mask_reori_scaled}"

        task_exec

        task_in="fslmaths ${fs_lobes_nii} -mas ${bmc_minL_conn} ${fs_lobes_minL_nii} && ImageMath 3 ${fs_lobes_plusL_nii} \
        addtozero ${fs_lobes_minL_nii} ${L_mask_reori_scaled}"

        task_exec

        task_in="labelconvert -force -nthreads ${ncpu} ${fs_parc_plusL_nii} ${function_path}/share/labelconvert/FreeSurferColorLUT+lesion.txt \
        ${function_path}/share/labelconvert/fs_default+lesion.txt ${fs_parc_plusL_nii_LC} && labelconvert -force -nthreads ${ncpu} \
        ${fs_parc_nii} ${FS_path1}/FreeSurferColorLUT.txt ${mrtrix_path}/share/mrtrix3/labelconvert/fs_default.txt ${fs_parc_nii_LC} \
        && labelconvert -force -nthreads ${ncpu} ${fs_parc_minL_nii} ${FS_path1}/FreeSurferColorLUT.txt \
        ${mrtrix_path}/share/mrtrix3/labelconvert/fs_default.txt ${fs_parc_minL_nii_LC}"

        task_exec
    
    else
        
        echo " lobes fs image already done, skipping. " | tee -a ${prep_log}
        
    fi

    # use for loop to read all values and indexes

    search_wf_mark6=($(find ${ROIs} -type f | grep LT_PWM_bin.nii.gz));
        
    if [[ -z "$search_wf_mark6" ]]; then

        for i in {0..11}; do

            echo "Now working on ${labels[$i]}" | tee -a ${prep_log}

            task_in="fslmaths ${fs_lobes_nii} -thr ${gm[$i]} -uthr ${gm[$i]} ${ROIs}/${labels[$i]}_gm.nii.gz"

            task_exec

            task_in="fslmaths ${fs_lobes_nii} -thr ${wm[$i]} -uthr ${wm[$i]} ${ROIs}/${labels[$i]}_wm.nii.gz"

            task_exec

            task_in="fslmaths ${ROIs}/${labels[$i]}_gm.nii.gz -add ${ROIs}/${labels[$i]}_wm.nii.gz -bin ${ROIs}/${labels[$i]}_bin.nii.gz"

            task_exec

        done

        i=""

        for i in {12..25}; do

            echo "Now working on ${labels[$i]}" | tee -a ${prep_log}

            task_in="fslmaths ${fs_lobes_nii} -thr ${gm[$i]} -uthr ${gm[$i]} -bin ${ROIs}/${labels[$i]}_bin.nii.gz"

            task_exec

        done

        i=""

        for i in {26..27}; do

            echo "Now working on ${labels[$i]}" | tee -a ${prep_log}

            task_in="fslmaths ${fs_lobes_nii} -thr ${wm[$i]} -uthr ${wm[$i]} -bin ${ROIs}/${labels[$i]}_bin.nii.gz"

            task_exec

        done
        
    else
        
        echo " isolating lobe labels already done, skipping to lesion overlap check" | tee -a ${prep_log}
        
    fi

    i=""

    # Now to check overlap and quantify existing overlaps
    # we also need to calculate volume and no. of vox for each lobe out of FS
    # also lesion volume

    l_vol=($(fslstats ${Lmask_o} -V))

    echo " * The lesion occupies " ${l_vol[0]} " voxels in total with " ${l_vol[0]} " cmm volume. " | tee -a $lesion_lobes_report

    for (( i=0; i<${labelslength}; i++ )); do


        task_in="fslmaths ${ROIs}/${labels[$i]}_bin.nii.gz -mas ${Lmask_o} ${overlap}/${labels[$i]}_intersect_L_mask.nii.gz"

        task_exec

        b=($(fslstats ${overlap}/${labels[$i]}_intersect_L_mask.nii.gz -V))
        
        a=($( echo ${b[0]} | cut -c1-1))

        vol_lobe=($(fslstats ${ROIs}/${labels[$i]}_bin.nii.gz -V))

        echo " - The " ${labels[$i]} " label is " ${vol_lobe[0]} " voxels in total, with a volume of " ${vol_lobe[1]} " cmm volume. " | tee -a ${lesion_lobes_report}

        if [[ $a -ne 0 ]]; then

            vol_ov=($(fslstats ${overlap}/${labels[$i]}_intersect_L_mask.nii.gz -V))
            
            ov_perc=($(echo "scale=4; (${vol_ov[1]}/${vol_lobe[1]})*100" | bc ))

            echo " ** The lesion overlaps with the " ${labels[$i]} " in " ${vol_ov[1]} \
            " cmm " ${ov_perc} " percent of total lobe volume " | tee -a ${lesion_lobes_report}

        else

        echo " No overlap between the lesion and " ${labels[$i]} " lobe. " | tee -a ${lesion_lobes_report}

        fi


    done

elif [[ "${P_flag}" -eq 0 ]] ; then

    echo
    echo "Fresurfer flag not set, finished, exiting" >&2
    echo "Fresurfer flag not set, finished, exiting" | tee -a ${prep_log}
    echo

fi


finish_t=$(date +%s)

# echo ${start_t}
# echo ${finish_t}

run_time_s=($(echo "scale=4; (${finish_t}-${start_t})" | bc ))
run_time_m=($(echo "scale=4; (${run_time_s}/60)" | bc ))
run_time_h=($(echo "scale=4; (${run_time_m}/60)" | bc ))

echo " execution took ${run_time_m} minutes, or approximately ${run_time_h} hours. " | tee -a ${prep_log}

# if not running FS, but MSBP should use something like this:
# to run MSBP after a recon-all run is finished
# 
# docker run -it --rm -v $(pwd)/BIDS:/bids_dir \
# -v $(pwd)/BIDS/derivatives:/output_dir \
# -v /NI_apps/freesurfer/license.txt:/opt/freesurfer/license.txt \
# sebastientourbier/multiscalebrainparcellator:v1.1.1 /bids_dir /output_dir participant \
# --participant_label PT_002 --isotropic_resolution 1.0 --thalamic_nuclei \
# --brainstem_structures --skip_bids_validator --fs_number_of_cores 12 \
# --multiproc_number_of_cores 12 2>&1 >> $(pwd)/MSBP_trial_run.txt
