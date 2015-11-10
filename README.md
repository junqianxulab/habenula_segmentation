# Habenula Segmentation
---
Semi-automated Human Habenula Segmentation Program

---
## Requirements
* Python: version 2.7
* Python libraries: numpy, scipy

---
## Input Files
* Currently, only nifti1 is the only supported file format.
* AC-PC aligned T1w image
* T2w image that is registered to T1w image
* Myelin map: T1w/T2w which can be generated using fslmaths in FSL: `fslmaths T1w -div T2w Myelin`
* Habenula center positions: a nifti1 file containing two markers - a voxel in the right habenula as 1 and a vxoel in the left habenula as 2
* Note: The PreFreeSurfer in Human Connectome Project (HCP) Pipelines (https://github.com/Washington-University/Pipelines) is  recommended for the AC-PC alignment and T2w-to-T1w registration.

---
## How to Run
* After reading all functions in segment_hb.py, call the segment_hb() function.

THIS PAGE IS UNDER CONSTRUCTION
