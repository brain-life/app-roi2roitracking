[![Abcdspec-compliant](https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg)](https://github.com/brain-life/abcd-spec)
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-brainlife.app.279-blue.svg)](https://doi.org/10.25663/brainlife.app.279)

# ROI to ROI Tracking 

This app will perform ensemble tracking between 2 or more cortical regions of interest (ROIs) from either a freesurfer parcellation or an atlas parcellation. This app uses MrTrix0.2.12's streamtrack probabilistic tracking functionality to track between pairs of ROIs. This app takes in the following required inputs: DWI, rois, Freesurfer. Optionally, a user can input precomputed CSDs, brainmask, and white matter masks. If these inputs are not provided, the app will generate these datatypes internally and output the derivatives. This is why the Freesurfer is a required input. If the CSD is computed internally, the derivatives will be outputted as a CSD datatype. If the user inputs the CSD, the files are copied over as a CSD datatype output.

To specify which ROIs are desired, the user can specify each pair as the ROI numbers found in the ROI datatype seperated by a space. For example, if the user wanted to track between ROIs 0001 and 0002 from their ROIs datatype, the user should input 0001 0002 in the roiPair field. The first value will generally be treated as the seed ROI and the second as the termination ROI. However, the user can specify seeding in both ROIs by setting the 'multiple_seeds' field to true. If the user wants to make multiple tracks, enter the user can input multiple pairs by creating a new line in the roiPair field. The output classification structure will then contain the same number of tracks as ROI pairs.

This app provides the user with a large number of exposed parameters to guide and shape tractography. These include a maximum spherical harmonic order (lmax), number of repetitions, minimum and maximum length of streamlines, step size, maximum number of attempts, streamline count, FOD cutoff, and minimum radius of curvature. For lmax, the user can specify whether or not to track on a single lmax or 'ensemble' across lmax's. If the user wants to track in just a single lmax, set the 'single_lmax' field to true. Else, leave as false. If the user does not know which lmax to use, they can leave the 'max_lmax' field empty and the app will automatically calculate the maximum lmax based on the number of unique directions in the non-b0 weighted volumes of the DWI. For FOD cutoff, minimum radius of curvature, and step size, the user can input multiple values to perform 'ensemble tracking'. If this is desired, the user can input each value separated by a space in the respective fields (example: 0.25 0.5). The outputs of each iteration will be merged together in the final output.

This app requires the ROIs and DWI datatypes to have the same dimensions. If the ROIs were generated with the 'Generate ROIs in dMRI space' app, then the ROIs and DWI will have the same dimensions. If another app was used to generate the ROIs input, the user will need to set the 'reslice' field to true. This will reslice the ROIs to have the same dimensions as the DWI image. This also assumes proper alignment between the DWI image and the ROIs.
This app uses multiple software packages, including Freesurfer, MrTrix0.2.12, Matlab, and python. 

### Authors 

- Brad Caron (bacaron@iu.edu)
- Ilaria Sani (isani01@rockefeller.edu) 

### Contributors 

- Soichi Hayashi (hayashis@iu.edu) 

### Funding 

[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-1636893](https://img.shields.io/badge/NSF_BCS-1636893-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-ACI-1916518](https://img.shields.io/badge/NSF_ACI-1916518-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1916518)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NIH-NIBIB-R01EB029272](https://img.shields.io/badge/NIH_NIBIB-R01EB029272-green.svg)](https://grantome.com/grant/NIH/R01-EB029272-01)

### Citations 

Please cite the following articles when publishing papers that used data, code or other resources created by the brainlife.io community. 

1. Fibre-tracking was performed using the MRtrix package (J-D Tournier, Brain Research Institute, Melbourne, Australia, http://www.brain.org.au/software/) (Tournier et al. 2012)
2. Tournier JD, Calamante F, Gadian DG, Connelly A Direct estimation of the fiber orientation density function from diffusion-weighted MRI data using spherical deconvolution Neuroimage 2004; 23 (3): 1176-1185
3. Tournier JD, Calamante F, Connelly A MRtrix: diffusion tractography in crossing fibre regions International Journal of Imaging Systems and Technology 2012; in press, DOI: 10.1002/ima.22005 

## Running the App 

### On Brainlife.io 

You can submit this App online at [https://doi.org/10.25663/brainlife.app.279](https://doi.org/10.25663/brainlife.app.279) via the 'Execute' tab. 

### Running Locally (on your machine) 

1. git clone this repo 

2. Inside the cloned directory, create `config.json` with something like the following content with paths to your input files. 

```json 
{
   "dwi":    "testdata/dwi/dwi.nii.gz",
   "bvals":    "testdata/dwi/dwi.bvals",
   "bvecs":    "testdata/dwi/dwi.bvals",
   "lmax2":    "testdata/csd/lmax2.nii.gz/",
   "rois":    "testdata/rois/rois/",
   "anat":    "testdata/anat/t1.nii.gz",
   "minlength":    10,
   "maxlength":    200,
   "max_lmax":    8,
   "cutoff":    "0.025",
   "roiPair":    "001 0002",
    "stepsize":    "0.25",
   "num_fibers":    500,
   "curv":    "45",
   "max_num":    100000,
   "single_lmax":    true,
   "reslice":    false,
   "multiple_seed":    false,
   "num_repetitions":    1
} 
``` 

### Sample Datasets 

You can download sample datasets from Brainlife using [Brainlife CLI](https://github.com/brain-life/cli). 

```
npm install -g brainlife 
bl login 
mkdir input 
bl dataset download 
``` 

3. Launch the App by executing 'main' 

```bash 
./main 
``` 

## Output 

The main output of this App is contains the whole-brain tractogram (tck) and the internally computed masks, a white-matter classification structure (wmc), and all of the other derivatives generated. If masks were inputted, the output is simply copies of the inputs. 

#### Product.json 

The secondary output of this app is `product.json`. This file allows web interfaces, DB and API calls on the results of the processing. 

### Dependencies 

This App requires the following libraries when run locally. 

- MrTrix0.2: 
- jsonlab: 
- singularity: 
- FSL: 
- Matlab 
