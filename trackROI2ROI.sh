#!/bin/bash

## Ilaria Sani - 20180313
## Human Engogenous Attention Network
## (modified after Ilaria Sani - 20170117 - NHP Endogenous Attention Network)
## (modified after Brent McPherson)
##
## MRTRIX PRE-PROCESSING (no need for WM mask for humans)
## This script process DTI data with MRTRIX, i.e. prepares the data for tracking
## it also creates a WM mask
## 
## NEEDS:
## 1. $dwi.mif
## 2. grads.b 
## 3. mask.mif 

########### DEFINE PATHS #####################
export PATH=$PATH:/usr/lib/mrtrix/bin

BGRAD="grad.b"
dtiinit=`jq -r '.dtiinit' config.json`
export input_nii_gz=$dtiinit/`jq -r '.files.alignedDwRaw' $dtiinit/dt6.json`
export BVALS=$dtiinit/`jq -r '.files.alignedDwBvals' $dtiinit/dt6.json`
export BVECS=$dtiinit/`jq -r '.files.alignedDwBvecs' $dtiinit/dt6.json`
fsurfer=`jq -r '.freesurfer' config.json`
rois=`jq -r '.parcellation' config.json`
roipair=`jq -r '.roiPair' config.json`
NUM=`jq -r '.num_fibers' config.json`
MAXNUM=`jq -r '.max_num' config.json`
STEPSIZE=`jq -r '.stepsize' config.json`
MINLENGTH=`jq -r '.minlength' config.json`
MAXLENGTH=`jq -r '.maxlength' config.json`
WMMK=wm_mask.mif


########### CREATE grad.b #####################
if [ -f grad.b ]; then
	echo "grad.b exist... skipping"
else
	echo "starting matlab to create grad.b"
    	./compiled/main
fi

#if max_lmax is empty, auto calculate
MAXLMAX=`jq -r '.max_lmax' config.json`
if [[ $MAXLMAX == "null" || -z $MAXLMAX ]]; then
    echo "max_lmax is empty... determining which lmax to use from .bvals"
    MAXLMAX=`./calculatelmax.py`
fi

########### CONVERSION OF FILES ############
mrconvert $input_nii_gz dwi.mif
mrconvert mask_anat.nii.gz b0.mif
mrconvert wm_anat.nii.gz $WMMK

mkdir roi;
for ROI in ${roipair[*]}
	do
		cp $rois/roi_${ROI}.nii.gz ./
		# add line to remove .nii.gz from name
		mrconvert roi_${ROI}.nii.gz roi_${ROI}.mif
		mv roi_${ROI}.nii.gz ./roi/
	done
	ret=$?	
	if [ ! $ret -eq 0 ]; then
		exit $ret
	fi

########### CREATE FILES FOR TRACKING ######
## create a t2-mask from b0
if [ -f mask.mif ]; then
	echo "t2-mask from b0 already exists...skipping"
else
	time average dwi.mif -axis 3 - | threshold - - | median3D - - | median3D - mask.mif
	ret=$?
	if [ ! $ret -eq 0 ]; then
		exit $ret
	fi
fi

## fit diffusion model
if [ -f dt.mif ]; then
	echo "diffusion tensor already exists...skipping"
else
	time dwi2tensor dwi.mif -grad $BGRAD dt.mif
	ret=$?
	if [ ! $ret -eq 0 ]; then
		exit $ret
	fi
fi

## create FA image
if [ -f fa.mif ]; then
	echo "FA image already exists...skipping"
else
	time tensor2FA dt.mif - | mrmult - mask.mif fa.mif
	ret=$?
	if [ ! $ret -eq 0 ]; then
		exit $ret
	fi
fi

## create rgb eingenvectors
if [ -f ev.mif ]; then
	echo "RGB eigenvectors already exists...skipping"
else
	time tensor2vector dt.mif - | mrmult - mask.mif ev.mif
	if [ ! $ret -eq 0 ]; then
		exit $ret
	fi
fi

## create single fiber mask
if [ -f sf.mif ]; then
	echo "Single fiber mask already exists...skipping"
else
	time erode mask.mif -npass 3 - | mrmult fa.mif - - | threshold - -abs 0.7 sf.mif
	if [ ! $ret -eq 0 ]; then
		exit $ret
	fi
fi

## create response numbers for CSD fit
for (( i_lmax=2; i_lmax<=$MAXLMAX; i_lmax+=2 )); do
	if [ -f response${i_lmax}.txt ]; then
    		echo "response${i_lmax}.txt already exist... skipping"
	else
    	time estimate_response -quiet dwi.mif sf.mif -grad $BGRAD -lmax $i_lmax response${i_lmax}.txt
    	ret=$?
    	if [ ! $ret -eq 0 ]; then
        exit $ret
    	fi
	fi
done

## fit csd model
for (( i_lmax=2; i_lmax<=$MAXLMAX; i_lmax+=2 )); do
# Perform CSD in each white matter voxel
	lmaxout=csd${i_lmax}.mif
	if [ -s $lmaxout ]; then
		echo "$lmaxout already exist - skipping csdeconv"
	else
		time csdeconv -quiet dwi.mif -grad $BGRAD response${i_lmax}.txt -lmax $i_lmax -mask mask.mif $lmaxout
		ret=$?
		if [ ! $ret -eq 0 ]; then
			exit $ret
		fi
	fi
done

################# ROI2ROI TRACKING ############################
ROI=(*roi_*.mif);
range=` expr ${#ROI[@]}`
nTracts=` expr ${range} / 2`
for (( i=0; i<=$nTracts; i+=2 ));
	do
		for i_track in 01 #02 03 04 05 06 07 08 09 10
			do
				echo ${i_track}
				for (( i_lmax=2; i_lmax<=$MAXLMAX; i_lmax+=2 ))
					do
    						for curv in 0.5 1 2 3 4
							do
								streamtrack SD_PROB csd${i_lmax}.mif tract$((i/2+1))_lmax${i_lmax}_crv${curv}_${i_track}.tck -mask $WMMK -grad $BGRAD -number $NUM -maxnum $MAXNUM -curvature $curv -step $STEPSIZE -minlength $MINLENGTH -length $MAXLENGTH -seed ${ROI[$((i))]} -include ${ROI[$((i))]} -include ${ROI[$((i+1))]} -stop
							done
					done
			done
		## concatenate tracts
		holder=(*tract$((i/2+1))*.tck)
		cat_tracks track$((i/2+1)).tck ${holder[*]}
		if [ ! $ret -eq 0 ]; then
			exit $ret
		fi
		rm -rf ${holder[*]}
		## tract info
		track_info track$((i/2+1)).tck > track_info$((i/2+1)).txt
		if [[ $((i/2+1)) == 1 ]];then
			mv track_info$((i/2+1)).txt track_info.txt
			mv track$((i/2+1)).tck track.tck
		fi
	done


################# CREATE CLASSIFICATION STRUCTURE ###############
./classification/classificationGenerator

################# CLEANUP #######################################
rm -rf ./roi/
rm -rf *.mif*
rm -rf grad.b
rm -rf *response*.txt
rm -rf *.nii.gz

