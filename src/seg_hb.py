#!/usr/bin/env python
# Habenula Segmentation
# 
# Arguments:
#   T1w
#   T2w
#
# Options:
#   Myelin
#   Hb_center
#   Warp_MNI
#   OutputBinary
#   OutputProb
#   Verbose
#   
#   MinVolume
#   T1Min
#   T1Max
#   T2Min
#   T2Max
#   Alpha
#   NumGrowing


import segment_hb
import argparse
import sys
import os
import subprocess

dirname = os.path.dirname(os.path.realpath(__file__))
hb_template_r = os.path.join(dirname, 'habenula_template_right.nii.gz')
hb_template_l = os.path.join(dirname, 'habenula_template_left.nii.gz')

def run_command(cmd, print_cmd=True):
    if print_cmd is True:
        print '>> %s' % cmd
    p = subprocess.call(cmd, shell=True)
    return p

def filename_wo_ext(filename):
    if filename[-7:] == '.nii.gz':
        return filename[:-7]

    elif filename[-4:] == '.nii':
        return filename[:-4]

    elif filename[-3:] == '.gz':
        return filename[:-3]

    return filename


parser = argparse.ArgumentParser()
parser.add_argument('-1', dest='t1w', help='T1w image', required=True)
parser.add_argument('-2', dest='t2w', help='T2w image', required=True)

parser.add_argument('-m', dest='myelin', default=None, help='Myelin image')

parser.add_argument('-c', dest='center', default=None, help='Habenula centers')
parser.add_argument('-w', dest='warp_to_acpc', default=None, help='Warp MNI to native')
parser.add_argument('-t', dest='template_threshold', default=0.2, type=float, help='Threshold for warped template')

parser.add_argument('-o', dest='output', default=None, help='Output binary filename')

parser.add_argument('-v', dest='verbose', type=int, default=0)
parser.add_argument('--min_volume', dest='min_volume', type=float, default=100)
parser.add_argument('--t1min', dest='t1min', default=None)
parser.add_argument('--t1max', dest='t1max', default=None)
parser.add_argument('--t2min', dest='t2min', default=10)
parser.add_argument('--t2max', dest='t2max', default=None)
parser.add_argument('--alpha', dest='alpha', type=float, default=0.9)
parser.add_argument('--num_growing', dest='growth', type=int, default=10)

ag = parser.parse_args()
if ag.center is None and ag.warp_to_acpc is None:
    sys.stderr.write('Either -c (habenula centers) or -w (warp to native) should be provided.\n')
    sys.exit(-1)

if ag.myelin is None:
    ag.myelin = filename_wo_ext(ag.t1w) + '_div_t2w.nii.gz'
    cmd = 'fslmaths %s -div %s %s' % (ag.t1w, ag.t2w, ag.myelin)
    run_command(cmd)

if ag.center is None:
    ag.center = filename_wo_ext(ag.myelin) + '_habenula_centers_from_template.nii.gz'
    tmp = 'tmp_hb_init.nii.gz'
    tmp_r = 'tmp_hb_init_r.nii.gz'
    tmp_l = 'tmp_hb_init_l.nii.gz'
    
    cmd = 'applywarp -i %s -r %s -o %s -w %s --interp=spline' % (hb_template_r, ag.t1w, tmp, ag.warp_to_acpc)
    run_command(cmd)
    cmd = 'fslmaths %s -thr %s -bin %s' % (tmp, ag.template_threshold, tmp_r)
    run_command(cmd)

    cmd = 'applywarp -i %s -r %s -o %s -w %s --interp=spline' % (hb_template_l, ag.t1w, tmp, ag.warp_to_acpc)
    run_command(cmd)
    cmd = 'fslmaths %s -thr %s -bin %s' % (tmp, ag.template_threshold, tmp_l)
    run_command(cmd)

    cmd = 'fslmaths %s -mul 2 -add %s %s' % (tmp_l, tmp_r, ag.center)
    run_command(cmd)

    cmd = 'rm -f %s %s %s' % (tmp, tmp_r, tmp_l)
    run_command(cmd)

if ag.output is None:
    ag.output = filename_wo_ext(ag.myelin) + '_habenula'

if ag.t1min is not None:
    ag.t1min = float(ag.t1min)
if ag.t1max is not None:
    ag.t1max = float(ag.t1max)
if ag.t2min is not None:
    ag.t2min = float(ag.t2min)
if ag.t2max is not None:
    ag.t2max = float(ag.t2max)

#print ag
#sys.exit(0)

return_values = segment_hb.segment_hb(
        ag.t1w,
        ag.t2w,
        ag.myelin,
        ag.center,
        ag.output,
        min_volume=ag.min_volume,
        t1min=ag.t1min,
        t1max=ag.t1max,
        t2min=ag.t2min,
        t2max=ag.t2max,
        sig_factor=ag.alpha,
        growth=ag.growth,
        verbose=ag.verbose
        )

print '%s %s' % return_values['volumes']

