#!/usr/bin/env python

import nibabel as nib
import numpy as np
import scipy.ndimage
import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--in', dest='fn_in', default=None, help='input probability segmentation filename', required=True)
parser.add_argument('--out', dest='bn_out', default=None, help='output basename, without .nii or .nii.gz', required=True)
parser.add_argument('--thr', dest='threshold', help='threshold, inclusive', type=float, default=0.5)
parser.add_argument('--each_hb', dest='each_hb', help='output each of left and right habenula', type=bool, default=True)

ag = parser.parse_args()

if not os.path.isfile(ag.fn_in):
    sys.stderr.write('ERROR: %s not found.\n' % ag.fn_in)
    sys.exit(1)

img = nib.load(ag.fn_in)
dat = img.get_fdata()

dat_bin = (dat >= ag.threshold).astype(np.int8)
label, n_label = scipy.ndimage.label(dat_bin, structure=np.ones((3,3,3)))

if n_label != 2:
    sys.stderr.write('ERROR: the number of connected objects are not 2 after applying threshold %s\n' % ag.threshold)
    sys.exit(2)

dat_l = (label==1).astype(np.int8)
dat_r = (label==2).astype(np.int8)

if (img.affine[0][0] > 0 and scipy.ndimage.center_of_mass(dat_l)[0] > scipy.ndimage.center_of_mass(dat_r)[0]) or \
   (img.affine[0][0] < 0 and scipy.ndimage.center_of_mass(dat_l)[0] < scipy.ndimage.center_of_mass(dat_r)[0]) :
    dat_l, dat_r = dat_r, dat_l

hdr_out = img.header.copy()
hdr_out.set_data_dtype(dat_l.dtype)

if ag.each_hb:
    for (lr, dat_seg) in [ ('_right', dat_r), ('_left', dat_l) ]:
        img_out = nib.Nifti1Image(dat_seg, img.affine, hdr_out)
        nib.save(img_out, ag.bn_out + lr + '.nii.gz')
else:
    img_out = nib.Nifti1Image(dat_l * 2 + dat_r, img.affine, hdr_out)
    nib.save(img_out, ag.bn_out + '.nii.gz')

