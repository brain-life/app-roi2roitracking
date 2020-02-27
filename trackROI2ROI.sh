#!/bin/bash

## Ilaria Sani - 20180313
## Human Engogenous Attention Network
## (modified after Ilaria Sani - 20170117 - NHP Endogenous Attention Network)
## (modified after Brent McPherson)
## (modified after Brad Caron)
## MRTRIX PRE-PROCESSING (no need for WM mask for humans)
## This script process DTI data with MRTRIX, i.e. prepares the data for tracking
## it also creates a WM mask
## 
## NEEDS:
## 1. $dwi.mif
## 2. grads.b 
## 3. mask.mif 

#make the script to fail if any of the command fails.
#set -e

##show commands executed (mainly for debugging)
set -x

########### DEFINE PATHS #####################
export PATH=$PATH:/usr/lib/mrtrix/bin

BGRAD="grad.b"
input_nii_gz=`jq -r '.dwi' config.json`
BVALS=`jq -r '.bvals' config.json`
BVECS=`jq -r '.bvecs' config.json`
brainmask=`jq -r '.brainmask' config.json`
dtiinit=`jq -r '.dtiinit' config.json`
fsurfer=`jq -r '.freesurfer' config.json`
rois=`jq -r '.rois' config.json`
roipair=`jq -r '.roiPair' config.json`
NUM=`jq -r '.num_fibers' config.json`
MAXNUM=`jq -r '.max_num' config.json`
STEPSIZE=`jq -r '.stepsize' config.json`
MINLENGTH=`jq -r '.minlength' config.json`
MAXLENGTH=`jq -r '.maxlength' config.json`
NUM_REPETITIONS=`jq -r '.num_repetitions' config.json`
CURVATURE=`jq -r '.curv' config.json`
MAXLMAX=`jq -r '.max_lmax' config.json`
lmax2=`jq -r '.lmax2' config.json`
lmax4=`jq -r '.lmax4' config.json`
lmax6=`jq -r '.lmax6' config.json`
lmax8=`jq -r '.lmax8' config.json`
lmax10=`jq -r '.lmax10' config.json`
lmax12=`jq -r '.lmax12' config.json`
lmax14=`jq -r '.lmax14' config.json`
single_lmax=`jq -r '.single_lmax' config.json`
multiple_seed=`jq -r '.multiple_seed' config.json`
WMMK=wm_mask.mif

# if dtiinit is inputted, set appropriate fields
if [[ ! ${dtiinit} == "null" ]]; then
	input_nii_gz=$dtiinit/*dwi_aligned*.nii.gz
	BVALS=$dtiinit/*.bvals
	BVECS=$dtiinit/*.bvecs
	brainmask=$dtiinit/dti/bin/brainMask.nii.gz
	mrconvert ${brainmask} mask.mif
fi

#generate grad.b from bvecs/bvals
#load bvals/bvecs
bvals=$(cat $BVALS)
bvecs_x=$(cat $BVECS | head -1)
bvecs_y=$(cat $BVECS | head -2 | tail -1)
bvecs_z=$(cat $BVECS | tail -1)
#convert strings to array of numbers
bvecs_x=($bvecs_x)
bvecs_y=($bvecs_y)
bvecs_z=($bvecs_z)
#output grad.b
i=0
true > grad.b
for bval in $bvals; do
    echo ${bvecs_x[$i]} ${bvecs_y[$i]} ${bvecs_z[$i]} $bval >> grad.b
    i=$((i+1))
done

#if max_lmax is empty, auto calculate

if [[ $MAXLMAX == "null" || -z $MAXLMAX ]]; then
    echo "max_lmax is empty... determining which lmax to use from .bvals"
    MAXLMAX=`./calculatelmax.py`
fi

if [ ! -f dwi.mif ]; then
    mrconvert $input_nii_gz dwi.mif
fi

if [ ! -f $WMMK ]; then
    mrconvert wm_anat.nii.gz $WMMK
fi

mkdir -p roi
re='^[0-9]+$'
for ROI in ${roipair[*]}
    do
        if ! [[ $ROI =~ $re ]]; then
	    if [ ! -f ./$ROI.nii.gz ]; then
            	cp $rois/$ROI.nii.gz ./
	    fi
	    if [ ! -f $roi_${ROI}.mif ]; then
                mrconvert $ROI.nii.gz roi_${ROI}.mif
            fi
            mv ${ROI}.nii.gz ./roi/
        else
	    if [ ! -f ./$ROI.nii.gz ]; then
            	cp $rois/ROI${ROI}.nii.gz ./
	    fi
            # add line to remove .nii.gz from name
            if [ ! -f roi_${ROI}.mif ]; then
                mrconvert ROI${ROI}.nii.gz roi_${ROI}.mif
            fi
            mv ROI${ROI}.nii.gz ./roi/
        fi
    done
	ret=$?	
	if [ ! $ret -eq 0 ]; then
		exit $ret
	fi

########### CREATE FILES FOR TRACKING ######
## create a t2-mask from b0
if [ -f mask.mif ]; then
	echo "t2-mask from b0 already exists...skipping"
elif [[ ${brainmask} == 'null' ]]; then
	time average dwi.mif -axis 3 - | threshold - - | median3D - - | median3D - mask.mif
	ret=$?
	if [ ! $ret -eq 0 ]; then
		exit $ret
	fi
else
	mrconvert ${brainmask} mask.mif
	ret=$?
	if [ ! $ret -eq 0 ]; then
		exit $ret
	fi
fi

# csd generation or copying
if [[ ${single_lmax} == true ]]; then
	lmaxs=$(seq ${MAXLMAX} ${MAXLMAX})
else
	lmaxs=$(seq 2 2 ${MAXLMAX})
fi

for LMAXS in ${lmaxs}; do
	# copy over lmax if inputted. else, create lmax
	if [ ! -f csd${MAXLMAX}.mif ]; then
		lmaxvar=$(eval "echo \$lmax${MAXLMAX}")
		if [[ ${lmaxvar} == 'null' ]]; then
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
			if [ -f response${LMAXS}.txt ]; then
			    	echo "response${LMAXS}.txt already exist... skipping"
			else
			    	time estimate_response -quiet dwi.mif sf.mif -grad $BGRAD -lmax ${LMAXS} response${LMAXS}.txt
			    	ret=$?
			    	if [ ! $ret -eq 0 ]; then
					exit $ret
			    	fi
			fi
			
			## fit csd model
			# Perform CSD in each white matter voxel
			lmaxout=csd${LMAXS}.mif
			if [ -s $lmaxout ]; then
				echo "$lmaxout already exist - skipping csdeconv"
			else
				time csdeconv -quiet dwi.mif -grad $BGRAD response${LMAXS}.txt -lmax $LMAXS -mask mask.mif $lmaxout
				ret=$?
				if [ ! $ret -eq 0 ]; then
					exit $ret
				fi
			fi
		else
			echo "csd already inputted. skipping csd generation"
			mrconvert ${lmaxvar} ./csd${MAXLMAX}.mif
		fi
	else
		echo "csd exists. skipping"
	fi
done

################# ROI2ROI TRACKING ############################
ROI=(*roi_*.mif);
range=` expr ${#ROI[@]}`
nTracts=` expr ${range} / 2`

if [[ ${multiple_seed} == true ]]; then
    mradd *roi_*.mif seed.mif
    seed=seed.mif
else
    seed_name="$(echo ${roipair} | cut -d' ' -f1)"
    seed=roi_${seed_name}.mif
fi

for (( i=0; i<=$nTracts; i+=2 )); do
	for LMAXS in ${lmaxs}; do
		for i_track in $(seq $NUM_REPETITIONS); do
			 echo ${i_track}
			 for curv in ${CURVATURE}; do
				 out=tract$((i/2+1))_lmax${LMAXS}_crv${curv}_${i_track}.tck
				 streamtrack -quiet SD_PROB csd${LMAXS}.mif tmp.tck \
				-mask $WMMK \
				-grad $BGRAD \
				-number $NUM \
				-maxnum $MAXNUM \
				-curvature $curv \
				-step $STEPSIZE \
				-minlength $MINLENGTH \
				-length $MAXLENGTH \
				-seed ${seed} \
				-include ${ROI[$((i))]} \
				-include ${ROI[$((i+1))]} \
				-stop
	    	    	        mv tmp.tck $out
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
./compiled/classificationGenerator

################# CLEANUP #######################################
if [ -f ./track/track/tck ]; then
	rm -rf ./roi/
	rm -rf *.mif*
	rm -rf grad.b
	rm -rf *response*.txt
else
	echo "tracking failed"
	exit 1
fi
#rm -rf *.nii.gz #this removes aparc+aseg.nii.gz and other .nii.gz needed later
