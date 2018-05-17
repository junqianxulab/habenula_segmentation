# Habenula Segmentation
---
Automated Human Habenula Segmentation Program

---
## Requirements
* Python: version 2.7
* Python libraries: numpy, scipy, nibabel
* FSL

---
## Input Files
* nifti1 is the only supported file format.
* AC-PC aligned T1w image
* T2w image that is registered to T1w image
* Habenula center positions:
    * a nifti1 file containing two markers - a voxel in the right habenula as 1 and a vxoel in the left habenula as 2
    * or, `MNI152_0.7mm-to-ACPC` warp file. (`MNINonLinear/xfms/standard2acpc_dc.nii.gz` if you used HCP pipeline)
* Optional: T1w/T2w which can be generated using fslmaths in FSL: `fslmaths T1w -div T2w T1wDivT2w`
* Note: The PreFreeSurfer in Human Connectome Project (HCP) Pipelines (https://github.com/Washington-University/Pipelines) is recommended for the AC-PC alignment and T2w-to-T1w registration.

---
## How to Run

### method 1: using command line

#### if you have habenula center position file, `center.nii.gz`

```
seg_hb.py -1 T1w.nii.gz -2 T2w.nii.gz  -c center.nii.gz [-m T1wDivT2w.nii.gz] [-o output.nii.gz]
```

#### if you have warp file, `standard2acpc_dc.nii.gz`

```
seg_hb.py -1 T1w.nii.gz -2 T2w.nii.gz -w standard2acpc_dc.nii.gz [-m T1wDivT2w.nii.gz] [-o output.nii.gz]
```

#### options

```
seg_hb.py -1 T1w.nii.gz -2 T2w.nii.gz [-m T1wDivT2w.nii.gz] [-o output.nii.gz] \
          [-c center.nii.gz] [-w warp_MNI152_0.7_to_ACPC] \
          [-t TEMPLATE_THRESHOLD] [-v VERBOSE] \
          [--min_volume MIN_VOLUME] [--t1min T1MIN] [--t1max T1MAX] \
          [--t2min T2MIN] [--t2max T2MAX] [--alpha ALPHA] \
          [--num_growing GROWTH]
```

---
### method 2: using python script
```
import segment
segment.segment_hb(
        'T1w_filename',
        'T2w_filename',
        'myelin_filename',
        'hb_center_filename',
        'output_filename',
        min_volume = 80,
        t1min = None,
        t1max = None,
        t2min = 10,
        t2max = None,
        sig_factor = 0.9,
        growth = 5,
        verbose = False,
        return_volume=True,
        return_step_volume=False,
        return_center_of_mass=False,
        return_hb_index=False,
        return_hb_neighbor_index=False
        )
```

#### input:
*  t1, t2, myelin: T1w, T2w, Myelin Nifti1 images
*  hb center: mask Nifti1 file indicating Hb centers (1: right Hb, 2: left Hb)

#### output:
*  output filename: output (segmentation) Nifti1 image file name

#### parameters:
*  `min_volume`: minimum ROI volume (mm^3) for each Hb. Default 80.
*  `t1min, t1max, t2min, t2max`: first thresold values for T1w, T2w.
*  `sig_factor`: Factor of sigma for myelin threshold for Hb segmentation. Default 0.9. The higher factor, the smaller Hb.
*  `growth`: Number of iteration for region growth method. 0 for disable it. Default 5.
*  `verbose`: If True, Display details. Default False.

#### return:
*  [color1 Hb volume, color2 Hb volume], [color1 Hb partial volume, color2 Hb partial volume], [list(left Hb index), list(right Hb index)]

# Habenula Template
We release habenula templates: `habenula_template` directory.
* `habenula_template_HCP_50_native_[left/right].nii.gz` : habenula template in the MNI space generated from HCP 50 subjects' habenula segmentation in the native space.

---
#### Reference:

[1] Kim et al, Human habenula segmentation using myelin content. Neuroimage, 2016, 130 : 145-156 http://www.ncbi.nlm.nih.gov/pubmed/26826517

[2] Kim et al, Reproducibility of myelin content‐based human habenula segmentation at 3 Tesla. Human Brain Mapping, 2018 https://doi.org/10.1002/hbm.24060

