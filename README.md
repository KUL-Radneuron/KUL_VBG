# KUL_VBG
A workflow to allow Freesurfer recon-all to run on brain images with large lesions.
VBG is a bash script tested in Mac OSX, Ubuntu 18.0 and CentOS. 

The first commit on this repository corresponds to the version of the workflow used in the preprint "Virtual brain grafting: Enabling whole brain parcellation in the presence of large lesions. Radwan et al., 2020, DOI: https://doi.org/10.1101/2020.09.30.20204701, available via: https://www.medrxiv.org/content/10.1101/2020.09.30.20204701v1). This work was published in Neuroimage, 2021 available here: https://doi.org/10.1016/j.neuroimage.2021.117731

**Updated Dependencies:**
a) ANTs v2.3.1 and ANTsX scripts
b) FSL v6.0
c) MRtrix3 v3.0.2-64-g3eadb340
d) HD-BET
f) Freesurfer v6.0

Inputs:

Obligatory: 
1- Input to -p flag (participant name in BIDS convetion, without the leading sub-). 
2- A nifti format T1 WI of a subject (input to -a flag)
3- Binary lesion mask (lesion = 1, background = 0) integer nifti format (input to -l flag)
4- Indicate lesion mask space (input to -z flag) N.B. the specified lesion mask must have the same dimensions and transform as the input T1 WI.

Optional:
1- Specify location of intermediate processing and output folders (-m and -o flags)
2- Specify number of parallel workers used (input to -n flag)
3- Specify type of filling (default = uVBG, to activate bVBG use the -t flag)
4- Specify age group of participant (default = adult, to activate pediatric friendly mode specify the -P flag)
5- To run Freesurfer recon-all after the lesion filling is finished, specify the -F flag.
6- Verbose mode = -v

Examples:

    - Using the unilateral VBG approach and HD-BET for brain extraction, input data is in BIDS format with only 1 session
    KUL_VBG.sh -p pat001 -b -n 6 -l /fullpath/lesion_T1w.nii.gz -z T1 -o /fullpath/output -B 1
    
    - Using the bilateral VBG approach and HD-BET for brain extraction, input data is not in BIDS, FreeSurfer is also called at the end
    KUL_VBG.sh -p pat001 -a /fullpath/sub-PT_T1w.nii.gz -n 6 -l /fullpath/lesion_T1w.nii.gz -z T1 -o /fullpath/output -t -B 1 -F
	

Purpose:

    The purpose of this workflow is to generate a lesion filled image, with healthy looking synthetic tissue in place of the lesion
    Essentially excising the lesion and grafting over the brain tissue defect in the MR image space
    

Required arguments:

    -p:  BIDS participant name (anonymised name of the subject without the "sub-" prefix)
    -b:  if data is in BIDS
    -l:  full path and file name to lesion mask file per session
    -z:  space of the lesion mask used (only T1 supported in this version)
    -a:  Input precontrast T1WIs


Optional arguments:

    -s:  session (of the participant)
    -t:  Use the VBG template to derive the fill patch (if used, template tissue is used alongside native tissue to create the donor brain)
    -E:  Treat as an extra-axial lesion (skip VBG bulk, fill lesion patch with 0s, run FS and subsequent steps)
    -B:  specify brain extraction method (1 = HD-BET, 2 = ANTs-BET), if not set ANTs-BET will be used by default
    -F:  Run Freesurfer recon-all, generate aparc+aseg + lesion and lesion report
    -P:  In case of pediatric patients - use pediatric template (NKI_under_10 in MNI)
    -m:  full path to intermediate output dir
    -o:  full path to output dir (if not set reverts to default output ./lesion_wf_output)
    -n:  number of cpu for parallelisation (default is 6)
    -v:  show output from mrtrix commands
    -h:  prints help menu

Notes: 

    - You can use -b and the script will find your BIDS files automatically
    - If your data is not in BIDS, then use -a without -b
    - This version is for validation only.
    - In case of trouble with HD-BET see lines 1124 - 1200)
    - cook_template_4VBG requires two brains with unilateral lesions on opposing sides
    - it is meant to facilitate the grafting process and minimize intensity differences
    - You need a high resolution T1 WI and a lesion mask in the same space for VBG to run
    - If you end up with an empty image, it is possible you have a mismatch between the T1 and lesion mask
