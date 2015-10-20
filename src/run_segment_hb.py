execfile('segment_hb.py')

list_filename = '../../hcp500/hblist'
out_dir = 'output/'
seg_filename = 'hb_segmentation'
out_file = 'hb_volumes.csv'

sig_factor = 0.9
min_volume = 100
growth = 10
t1min = None
t1max = None
t2min = 10
t2max = None
verbose = True
#verbose = False

fin = open(list_filename)
hblist = fin.readlines()
hblist=[ hb[:-1] for hb in hblist if hb[0] != '#' ]
fin.close()

fout = open(out_file, 'w')
fout.write('Subject,Threshold_R,Threshold_L,Region_R,Region_L,Geometric_R,Geometric_L,Partial_R,Partial_L\n')

mkdir(out_dir)
for filename in hblist:
    print ''
    print 'Subject = ', filename
    base = '../../hcp500/data/' + filename + '/T1w/'
    return_values = segment_hb(
            base+'T1w_acpc_dc.nii.gz',
            base+'T2w_acpc_dc.nii.gz',
            base+'Myelin.nii.gz',
            base+'hb_center.nii.gz',
            out_dir+filename+'_'+seg_filename,
            min_volume=min_volume,
            t1min=t1min, t1max=t1max,
            t2min=t2min, t2max=t2max,
            sig_factor=sig_factor,
            growth=growth,
            verbose=verbose,
            return_step_volume=True )

    fout.write('%s' % filename)
    for volume_lr in return_values['step_volumes']:
        fout.write(',%s,%s' % (volume_lr[0], volume_lr[1]))
    fout.write('\n')
fout.close()

