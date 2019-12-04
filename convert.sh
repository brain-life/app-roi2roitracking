#!/bin/bash

## Create white matter mask and move rois to diffusion space for tracking

#exit if any command fails
#set -e 

#show commands runnings
#set -x

input_nii_gz=$(jq -r .dwi config.json)
bvecs=`jq -r '.bvecs' config.json`
bvals=`jq -r '.bvals' config.json`
dtiinit=`jq -r '.dtiinit' config.json`
fsurfer=`jq -r '.freesurfer' config.json`
eccentricity=`jq -r '.eccentricity' config.json`
hemi="lh rh"

if [[ ! ${dtiinit} == "null" ]]; then
        export input_nii_gz=$dtiinit/`jq -r '.files.alignedDwRaw' $dtiinit/dt6.json`
fi

mri_label2vol --seg $fsurfer/mri/aparc+aseg.mgz \
    --temp $input_nii_gz \
    --regheader $fsurfer/mri/aparc+aseg.mgz \
    --o aparc+aseg.nii.gz
    
mri_binarize --i aparc+aseg.nii.gz --min 1 --o mask_anat.nii.gz 

mri_binarize --i aparc+aseg.nii.gz --o wm_anat.nii.gz \
	--match 2 41 16 17 28 60 51 53 12 52 13 18 54 50 11 251 252 253 254 255 10 49 46 7

for HEMI in $hemi
do
        mri_vol2vol --mov $fsurfer/mri/${HEMI}.ribbon.mgz --targ ${input_nii_gz} --regheader --o ${HEMI}.ribbon.nii.gz
done

mri_vol2vol --mov ${fsurfer}/mri/ribbon.mgz --targ ${input_nii_gz} --regheader --o ribbon.nii.gz

mri_vol2vol --mov ${fsurfer}/mri/aparc.a2009s+aseg.mgz --targ ${input_nii_gz} --regheader --o aparc.a2009s.aseg.nii.gz

mri_vol2vol --mov ${eccentricity} --targ ${input_nii_gz} --regheader --o eccentricity.nii.gz
