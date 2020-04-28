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


########### DEFINE PATHS #####################
export PATH=$PATH:/usr/lib/mrtrix/bin

BGRAD="grad.b"
#input_nii_gz=`jq -r '.dwi' config.json`
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
response=`jq -r '.response' config.json`
single_lmax=`jq -r '.single_lmax' config.json`
multiple_seed=`jq -r '.multiple_seed' config.json`
white_matter=`jq -r '.white_matter' config.json`
flip_lr=`jq -r '.flip_lr' config.json`
WMMK=wm_mask.mif

mkdir -p csd track 

set -x
set -e

# if dtiinit is inputted, set appropriate field
if [[ ! ${dtiinit} == "null" ]]; then
    brainmask=$dtiinit/dti/bin/brainMask.nii.gz
    mrconvert ${brainmask} mask.mif
fi

#if max_lmax is empty, auto calculate
if [[ $MAXLMAX == "null" || -z $MAXLMAX ]]; then
    echo "max_lmax is empty... determining which lmax to use from .bvals"
    MAXLMAX=`./calculatelmax.py`
fi

#if [ ! -f dwi.mif ]; then
#    mrconvert $input_nii_gz dwi.mif
#fi

if [ ! -f $WMMK ]; then
    if [[ ${white_matter} == 'null' ]]; then
        mrconvert wm_anat.nii.gz $WMMK 
    else
        mrconvert ${white_matter} $WMMK
    fi
fi

#convert requested rois to mif
re='^[0-9]+$'
for ROI in ${roipair[*]}; do
echo "converting $ROI to mif........"
    if [[ $ROI =~ $re ]]; then
        echo "ROI is all numbers (assume ROI<num>.nii.gz format)"
        #cp $rois/ROI${ROI}.nii.gz ./
        # add line to remove .nii.gz from name
        if [ ! -f roi_${ROI}.mif ]; then
            mrconvert $rois/ROI${ROI}.nii.gz roi_${ROI}.mif
        fi
        #mv ROI${ROI}.nii.gz ./roi/
    else
        echo "ROI is not all numbers (assume <name>.nii.gz)"
        #cp $rois/$ROI.nii.gz ./
        if [ ! -f roi_${ROI}.mif ]; then
            mrconvert $rois/$ROI.nii.gz roi_${ROI}.mif 
        fi
        #mv ${ROI}.nii.gz ./roi/
    fi
done

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
            if [ ! -f ./csd/lmax${LMAXS}.nii.gz ]; then
                mrconvert csd${LMAXS}.mif ./csd/lmax${LMAXS}.nii.gz
                ret=$?
                if [ ! $ret -eq -0 ]; then
                    exit $ret
                fi
            fi
            if [[ ${LMAXS} == ${MAXLMAX} ]]; then
                cp response${LMAXS}.txt ./csd/response.txt
            fi
        else
            echo "csd already inputted. skipping csd generation"
            mrconvert ${lmaxvar} ./csd${MAXLMAX}.mif
            cp -v ${lmaxvar} ./csd/
            cp -v ${response} ./csd/ 
        fi
    else
        echo "csd exists. skipping"
    fi
done

################# ROI2ROI TRACKING ############################
#ROI=(roi_*.mif);
pairs=($roipair)
range=` expr ${#pairs[@]}`
nTracts=` expr ${range} / 2`

for (( i=0; i<=$nTracts; i+=2 )); do
    echo "creating seed for tract $((i/2))"
    roi1=roi_${pairs[$((i))]}.mif
    roi2=roi_${pairs[$((i+1))]}.mif
    if [[ ${multiple_seed} == true ]]; then
        #combine 2 rois into a single seed
        seed=seed_${pairs[$((i))]}_${pairs[$((i+1))]}.mif
        [ ! -f $seed ] && mradd $roi1 $roi2 $seed
    else
        #pick the first one for seed
        seed=roi_${pairs[$((i))]}.mif
    fi

    for LMAXS in ${lmaxs}; do
        for i_track in $(seq $NUM_REPETITIONS); do
             for curv in ${CURVATURE}; do
                 echo "tract $((i/2+1)) of $((nTracts/2)) / LMAXS:$LMAXS / repetition:$i_track / curv:$curv ------------------------"
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
                    -include $roi1 \
                    -include $roi2 \
                    -stop
                mv tmp.tck $out
            done
        done
    done

    ## concatenate tracts
    holder=(*tract$((i/2+1))*.tck)
    cat_tracks track$((i/2+1)).tck ${holder[*]}
    rm -rf ${holder[*]}
    
    ## tract info
    track_info track$((i/2+1)).tck > track_info$((i/2+1)).txt
done

echo "done with tracking"

