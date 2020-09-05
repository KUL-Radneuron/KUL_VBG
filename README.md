# KUL_VBG
A workflow to allow Freesurfer recon-all to run on brain images with large lesions.
VBG is a bash script tested in Mac OSX, Ubuntu 18.0 and CentOS.

Dependencies:
a) ANTs and ANTsX scripts
b) FSL
c) MRtrix3
d) HD-BET
f) Freesurfer

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

