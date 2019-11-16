###########################################################################
# Habeluna Segmentation
# 
# Modified: Aug/07/2015
###########################################################################

# dependency: nibabel, numpy, scipy

import errno
import os
import numpy as np
from scipy.optimize import curve_fit
from scipy import ndimage

def mkdir(directory):
    '''make a directory
    '''
    try:
        os.makedirs(directory)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(directory):
            pass
        else: raise

class FindNbd:
    '''Class to find neighbors inside boundary'''
    __slots__ = ['nx', 'ny', 'nz']
    def __init__(self, shape):
        self.nx, self.ny, self.nz = shape
    def get_nbds_26(self, index):
        '''return 26-connected neighbors'''
        return [ (index[0]+i,index[1]+j,index[2]+k)    \
                        for i in [-1,0,1] for j in [-1,0,1] for k in [-1,0,1] \
                        if not i == j == k == 0 \
                        if 0 <= index[0]+i < self.nx \
                        if 0 <= index[1]+j < self.ny \
                        if 0 <= index[2]+k < self.nz \
                        ]
    def get_nbds_6(self, index):
        '''return 6-connected neighbors'''
        return [ (index[0]+i,index[1]+j,index[2]+k)    \
                        for i,j,k in \
                            [ (-1,0,0), (0,-1,0), (0,0,-1), \
                                (1,0,0), (0,1,0), (0,0,1) ] \
                        if 0 <= index[0]+i < self.nx \
                        if 0 <= index[1]+j < self.ny \
                        if 0 <= index[2]+k < self.nz \
                        ]

class Options:
    '''Class containing segmentation parameters and data
    init:
            data_t1, data_t2, data_my, : numpy ndarray
            seed_voxels: right and left habenula seed voxels. e.g [ [(1,2,3),(1,2,4)], [(5,6,7)] ]
            zooms: resolution

            out_partial=None
            out_binary=None
            min_volume=80: minimum volume for each habenula's roi
            t1min=None
            t1max=None
            t2min=0
            t2max=None
            sigma_factor=0.9
            region_growth_maxiter=5
            template=None
            verbose=False
    '''
    __slots__ = ['data_t1', 'data_t2', 'data_my', 'seed_voxels', 'out_partial', 'out_binary',
            'shape', 'zooms', 'unit', 'get_nbd',
            'min_volume', 't1min', 't1max', 't2min', 't2max', 'sigma_factor',
            'region_growth_maxiter', 'template', 'verbose']
    def __init__(self, data_t1, data_t2, data_my, seed_voxels, zooms,
            out_partial=None,
            out_binary=None,
            min_volume=80,
            t1min=None,
            t1max=None,
            t2min=0,
            t2max=None,
            sigma_factor=0.9,
            region_growth_maxiter=5,
            template=None,
            verbose=False
            ):
        self.data_t1 = data_t1
        self.data_t2 = data_t2
        self.data_my = data_my
        self.seed_voxels = seed_voxels
        self.shape = data_t1.shape

        self.zooms = zooms
        self.unit = self.zooms[0] * self.zooms[1] * self.zooms[2]
        self.get_nbd = FindNbd(self.shape)
    
        self.min_volume = min_volume
        self.t1min = t1min
        self.t1max = t1max
        self.t2min = t2min
        self.t2max = t2max
        self.sigma_factor = sigma_factor
        self.region_growth_maxiter = region_growth_maxiter
        self.template = template
        self.verbose = verbose

        if out_partial is not None:
            self.out_partial = out_partial
        else:
            self.out_partial = np.zeros( data_t1.shape, dtype='<f4' )
        if out_binary is not None:
            self.out_binary = out_binary
        else:
            self.out_binary = np.zeros( data_t1.shape, dtype='<i2' )

    def print_options(self):
        print 'Options'

def histogram_gaussian_fitting(data, bins, normed=False, name='', rangex=None, fit_min=None, fit_max=None, verbose=False, weight=''):
    y, xx = np.histogram(data, bins=bins, normed=normed, range=rangex)
    x = np.array([ (xx[i]+xx[i+1])/2.0 for i in range(len(xx)-1) ])
   
    if fit_min == min(data) or fit_max == max(data):
        argmax_y = np.argmax(y)
    else:
        if weight == 'T1':
            argmax_y = len(y)/2 + np.argmax(y[len(y)/2:])
        elif weight == 'T2':
            argmax_y = np.argmax(y[:len(y)/2])
        else:
            argmax_y = np.argmax(y)
    m = x[argmax_y]
    n = len(x)

    x_fit = x[:]
    y_fit = y[:]
    if fit_min=='auto' or (weight == 'T1' and fit_min is None):
        x_fit = x_fit[int((3*argmax_y-n+1)/2):]
        y_fit = y_fit[int((3*argmax_y-n+1)/2):]
    elif type(fit_min) == type(1) or type(fit_min) == type(1.0):
        x_fit = x_fit[x_fit>fit_min]
        y_fit = y_fit[x_fit>fit_min]
    if fit_max=='auto' or (weight == 'T2' and fit_max is None):
        x_fit = x_fit[:int(3*argmax_y/2)]
        y_fit = y_fit[:int(3*argmax_y/2)]
    elif type(fit_max) == type(1) or type(fit_max) == type(1.0):
        x_fit = x_fit[x_fit<fit_max]
        y_fit = y_fit[x_fit<fit_max]

    try:
        popt, pcov = curve_fit(gaus,x_fit[y_fit>0],y_fit[y_fit>0],p0=[np.max(y_fit),m,np.std(data)])
    except RuntimeError:
        a, avg, sig = np.max(y), m, np.std(data)
        if verbose: print 'error while fitting in %s' % name
    except:
        a, avg, sig = np.max(y), m, np.std(data)
        if verbose: print 'error while fitting in %s' % name
    else:
        a, avg, sig = popt

    sig = abs(sig)
    return xx, y, a, avg, sig

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def read_nifti1_files(t1, t2, mye):
    import nibabel as nib
    nif_1 = nib.load(t1)
    header = nif_1.get_header()
    affine = nif_1.get_affine()
    zooms = header.get_zooms()
    
    nif_2 = nib.load(t2)
    if nif_2.shape != nif_1.shape:
        print 'shape of T2 does not match'
        os.exit(-1)
    if nif_2.get_header().get_zooms() != zooms:
        print 'resolution of T2 does not match'
        os.exit(-1)

    nif_3 = nib.load(mye)
    if nif_3.shape != nif_1.shape:
        print 'shape of Myelin does not match'
        os.exit(-1)
    if nif_3.get_header().get_zooms() != zooms:
        print 'resolution of Myelin does not match'
        os.exit(-1)

    dat_1 = nif_1.get_data()
    dat_2 = nif_2.get_data()
    dat_3 = nif_3.get_data()
    # swap if affine[0][0] > 0
    if affine[0][0] > 0:
        dat_1_lr = np.empty(dat_1.shape, dtype=dat_1.dtype)
        dat_2_lr = np.empty(dat_1.shape, dtype=dat_2.dtype)
        dat_3_lr = np.empty(dat_1.shape, dtype=dat_3.dtype)
        for z in range(dat_1.shape[2]):
            dat_1_lr[:,:,z] = np.flipud(dat_1[:,:,z])
            dat_2_lr[:,:,z] = np.flipud(dat_2[:,:,z])
            dat_3_lr[:,:,z] = np.flipud(dat_3[:,:,z])
    else:
        dat_1_lr = dat_1
        dat_2_lr = dat_2
        dat_3_lr = dat_3

    return dat_1_lr, dat_2_lr, dat_3_lr, affine, header

def write_nifti1_file(out_dat, out_filename, affine, header):
    import nibabel as nib
    header.set_data_dtype(out_dat.dtype)

    if affine[0][0] > 0:
        out_dat_lr = np.empty(out_dat.shape, dtype=out_dat.dtype)
        for z in range(out_dat.shape[2]):
            out_dat_lr[:,:,z] = np.flipud(out_dat[:,:,z])
    else:
        out_dat_lr = out_dat
    out_img = nib.Nifti1Image(out_dat_lr, affine, header)
    nib.save(out_img, out_filename)

def set_min_max(avg1, sig1, avg2, sig2, pt1min, pt1max, pt2min, pt2max):
    if pt1min is None:
        t1min = max([0, avg1 - 2*sig1])
    elif pt1min == 'avg':
        t1min = avg1
    else:
        t1min = pt1min
    if pt1max is None:
        t1max = np.inf
    else:
        t1max = pt1max

    if pt2max is None:
        #t2max = xx[-1] + 1
        t2max = avg2 + 2*sig2
    elif pt2max == 'avg':
        t2max = avg2
    else:
        t2max = pt2max

    return t1min, t1max, pt2min, t2max

def dilate_volume(base, get_nbd, volume, factor):
    '''Add nbds until >= volume'''
    done = base[:]
    while len(done)*factor < volume:
        to_do = []
        for ijk in done:
            nbds = get_nbd.get_nbds_6(ijk)
            for nbd in nbds:
                if nbd in done:
                    continue
                if nbd in to_do:
                    continue
                to_do.append(nbd)
        done = done + to_do
    return done

def dilate_intensity(base, domain, dat, get_nbd, intensity, conn=6):
    '''Add nbds >= intensity'''
    hb = base[:]
    while(True):
        to_do = []
        for ijk in hb:
            if conn == 26:
                nbds = get_nbd.get_nbds_26(ijk)
            else:
                nbds = get_nbd.get_nbds_6(ijk)

            for nbd in nbds:
                if nbd not in domain:
                    continue
                if nbd in to_do:
                    continue
                if nbd in hb:
                    continue
                if dat[nbd] < intensity:
                    continue

                to_do.append(nbd)

        if to_do == []:
            break
        hb = hb + to_do
    return hb

def dilate_out(hb, domain, get_nbd):
    '''dilate once within domain'''
    hb_out = hb[:]
    to_do = []
    for ijk in hb_out:
        nbds = get_nbd.get_nbds_6(ijk)
        for nbd in nbds:
            if nbd not in domain:
                continue
            if nbd in to_do:
                continue
            if nbd in hb_out:
                continue

            to_do.append(nbd)
            
    hb_out = hb_out + to_do
    return hb_out

def dilate_once(hb, get_nbd):
    '''dilate once'''
    hb_out = hb[:]
    to_do = []
    for ijk in hb_out:
        nbds = get_nbd.get_nbds_6(ijk)
        for nbd in nbds:
            if nbd in to_do:
                continue
            if nbd in hb_out:
                continue

            to_do.append(nbd)
            
    hb_out = hb_out + to_do
    return hb_out

def set_bdry(inside, outside, get_nbd):
    bdry = []
    for ijk in inside:
        nbds = get_nbd.get_nbds_6(ijk)
        for nbd in nbds:
            if (outside is None or nbd in outside) and ijk not in bdry:
                bdry.append(ijk)
    return bdry

def set_bdry_out(inside, outside, get_nbd):
    bdry = []
    for ijk in inside:
        nbds = get_nbd.get_nbds_6(ijk)
        for nbd in nbds:
            if (outside is None or nbd in outside) and nbd not in inside and nbd not in bdry:
                bdry.append(nbd)
    return bdry

def calculate_pivot(hb, hb_out, dat, factor = 1.0):
    std_factor = 3.0
    int_in = [ dat[ijk] for ijk in hb ]
    int_out = [ dat[ijk] for ijk in hb_out ]

    if int_out == []:
        return 0

    in_mean = np.mean(int_in)
    in_std = np.std(int_in)
    #in_mean = np.mean([value for value in int_in if in_mean-std_factor*in_std<value<in_mean+std_factor*in_std])
    in_mean = np.mean([value for value in int_in if value<in_mean+std_factor*in_std])
    #in_std = np.std([value for value in int_in if in_mean-std_factor*in_std<value<in_mean+std_factor*in_std])
    in_std = np.std([value for value in int_in if value<in_mean+std_factor*in_std])
    out_mean = np.mean(int_out)
    out_std = np.std(int_out)
    #out_mean = np.mean([value for value in int_out if out_mean-std_factor*out_std<value<out_mean+std_factor*out_std])
    #out_std = np.std([value for value in int_out if out_mean-std_factor*out_std<value<out_mean+std_factor*out_std])
    if in_mean < out_mean:
        print 'error: intensity of Hb < intensity of bdry'
        return 0
    if in_std == 0:
        print 'warning: std(Hb) = 0'
        return 0
    if out_std == 0:
        print 'warning: std(Hb_out) = 0'
        return 0
    in_istd = np.std(int_in)
    #in_istd = 1.0/in_istd
    #FIXME
    in_istd = 1.0/in_std
    out_istd = 1.0/out_std

    # factor*(pivot - out_mean) / out_std = (in_mean - pivot) / in_std
    pivot = (in_mean*in_istd + factor*out_mean*out_istd)/(in_istd + factor*out_istd)
    #FIXME
    #pivot = min([pivot, out_mean+out_std])

    return pivot

def refine_hb(hb, hb_out, dat, get_nbd, up_limit=0, verbose=False):
    stat_prev_mean_i = np.mean([dat[ijk] for ijk in hb])
    stat_prev_mean_o = np.mean([dat[ijk] for ijk in hb_out])

    bdry_in = set_bdry(hb, hb_out, get_nbd)
    bdry_out = set_bdry(hb_out, hb, get_nbd)

    pivot_up = calculate_pivot(hb, hb_out, dat, 1.0)
    pivot_down = calculate_pivot(hb, hb_out, dat, 1.0)
    #pivot_up = calculate_pivot(hb, hb_out, dat, 1.0/1.5)
    #pivot_down = calculate_pivot(hb, hb_out, dat, 1.5)
    if pivot_up == 0 or pivot_down == 0:
        return hb[:], (0,0)

    if pivot_up < up_limit:
        pivot_up = up_limit

    if pivot_up < pivot_down:
        return hb[:], (0,0)

    #pivot_up = calculate_pivot(hb, hb_out, dat, 1.0)
    #pivot_down = calculate_pivot(hb, hb_out, dat, 1.0)
    if verbose: print '  pivot up=%s, down=%s' % (pivot_up, pivot_down)

    # add
    hb_re = hb[:]
    num_up = 0
    for ijk in bdry_out:
        if dat[ijk] > pivot_up:
            if ijk not in hb_re:
                hb_re.append(ijk)
                num_up += 1

    hb_out_re = [ijk for ijk in hb_out if ijk not in hb_re]
    
    bdry_in = set_bdry(hb_re, hb_out_re, get_nbd)
    bdry_out = set_bdry(hb_out_re, hb_re, get_nbd)

    # remove
    num_down = 0
    for ijk in bdry_in:
        if dat[ijk] < pivot_down:
            if ijk not in hb_out_re:
                hb_out_re.append(ijk)
                num_down += 1

    hb_re = [ijk for ijk in hb_re if ijk not in hb_out_re]
    
    stat_next_mean_i = np.mean([dat[ijk] for ijk in hb_re])
    stat_next_mean_o = np.mean([dat[ijk] for ijk in hb_out_re])

    #if verbose: print '  Hb refine: in %s, out %s, mean_in %.2f -> %.2f, mean_out %.2f -> %.2f, mean_diff %.2f -> %.2f' \
    #        % (num_up, num_down, stat_prev_mean_i, stat_next_mean_i, stat_prev_mean_o, stat_next_mean_o, stat_prev_mean_i-stat_prev_mean_o, stat_next_mean_i-stat_next_mean_o)
    return hb_re, (num_up, num_down)

def ijk_add(ijk, list_add):
    return [ (ijk[0] + ijk_add[0], ijk[1] + ijk_add[1], ijk[2] + ijk_add[2]) for ijk_add in list_add ]

def get_a_partial_volume(ijk, hb, bdry_in, bdry_out, dat, ratio_nbd, get_nbd, verbose, use_ratio=False):
    int_in = 0.0
    int_out = 0.0
    num_in = 0.0
    num_out = 0
    factor_bdry = 0.5
    nbd_1 = [ [1,0,0], [0,1,0], [0,0,1], [-1,0,0], [0,-1,0], [0,0,-1] ]
    nbd_sqrt2 = [ [0,1,1], [0,1,-1], [0,-1,1], [0,-1,-1], [1,0,1], [1,0,-1], [1,1,0], [1,-1,0], [-1,0,1], [-1,0,-1], [-1,1,0], [-1,-1,0] ]
    nbd_sqrt3 = [ [1,1,1], [1,1,-1], [1,-1,1], [1,-1,-1], [-1,1,1], [-1,1,-1], [-1,-1,1], [-1,-1,-1] ]
    nbd_2 = [ [2,0,0], [0,2,0], [0,0,2], [-2,0,0], [0,-2,0], [0,0,-2] ]

    #if True:
    if False:
        for pairs in [ (nbd_1, 1.0), (nbd_sqrt2, 1.0/np.sqrt(2)), (nbd_sqrt3, 1.0/np.sqrt(3)), (nbd_2, 0.5) ]:
            nbds = ijk_add(ijk, pairs[0])

            for nbd in nbds:
                if nbd in bdry_in:
                    int_in += dat[nbd] * pairs[1] * factor_bdry
                    num_in += pairs[1] * factor_bdry
                if nbd in hb:
                    int_in += dat[nbd] * pairs[1]
                    num_in += pairs[1]
                elif nbd in bdry_out:
                    int_out += dat[nbd] * pairs[1] * factor_bdry
                    num_out += pairs[1] * factor_bdry
                else:
                    int_out += dat[nbd] * pairs[1]
                    num_out += pairs[1]

    # basic method
    #if True:
    if False:
        #nbds = get_nbd.get_nbds_26(ijk)
        nbds = ijk_add(ijk, nbd_1 + nbd_sqrt2 + nbd_sqrt3 + nbd_2)
        for nbd in nbds:
            if nbd in bdry_in:
                int_in += dat[nbd] * factor_bdry
                num_in += factor_bdry
            elif nbd in hb:
                int_in += dat[nbd]
                num_in += 1
            elif nbd in bdry_out:
                int_out += dat[nbd] * factor_bdry
                num_out += factor_bdry
            else:
                int_out += dat[nbd]
                num_out += 1

    #if True:
    if False:
        nbds = get_nbd.get_nbds_26(ijk)
        #nbds = ijk_add(ijk, nbd_1 + nbd_sqrt2 + nbd_sqrt3 + nbd_2)
        for nbd in nbds:
            if nbd in bdry_in:
                int_in += dat[nbd] * ratio_nbd[nbd]
                num_in += ratio_nbd[nbd]
            elif nbd in hb:
                int_in += dat[nbd]
                num_in += 1
            elif nbd in bdry_out:
                int_out += dat[nbd] * ratio_nbd[nbd]
                num_out += ratio_nbd[nbd]
            else:
                int_out += dat[nbd]
                num_out += 1

    if True:
    #if False:
        nbds = get_nbd.get_nbds_26(ijk)
        #nbds = ijk_add(ijk, nbd_1 + nbd_sqrt2 + nbd_sqrt3 + nbd_2)
        for nbd in nbds:
            if nbd in bdry_in:
                if use_ratio:
                    int_in += ratio_nbd[nbd]
                    num_in += 1
                else:
                    int_in += dat[nbd] * factor_bdry
                    num_in += factor_bdry
            elif nbd in hb:
                int_in += dat[nbd]
                num_in += 1
            elif nbd in bdry_out:
                if use_ratio:
                    int_out += ratio_nbd[nbd]
                    num_out += 1
                else:
                    int_out += dat[nbd] * factor_bdry
                    num_out += factor_bdry
            else:
                int_out += dat[nbd]
                num_out += 1

    if num_in == 0:
        return 1.0
    if num_out == 0:
        return 1.0

    int_in  /= num_in
    int_out /= num_out
    
    value = dat[ijk]
    if value >= int_in:
        return 1.0
    if value <= int_out:
        return 0.0
    if int_in < int_out:
        if verbose: print 'i(in) < i(out) at ', ijk
        return 0.5

    if not use_ratio:
        ratio_nbd[ijk] = int_in

    partial_volume = (value - int_out) / (int_in - int_out)
    if partial_volume > 1:
        partial_volume = 1.0
    elif partial_volume < 0:
        partial_volume = 0.0

    return partial_volume

def count_nbds(hb, bdry_in, bdry_out, get_nbd):
    ratio_nbd = {}
    for ijk in bdry_in:
        num_in = 0
        num_out = 0

        nbds = get_nbd.get_nbds_26(ijk)
        for nbd in nbds:
            if nbd in hb:
                num_in += 1
            else:
                num_out += 1
        ratio_nbd[ijk] = float(num_in)/(num_in+num_out)

    for ijk in bdry_out:
        num_in = 0
        num_out = 0

        nbds = get_nbd.get_nbds_26(ijk)
        for nbd in nbds:
            if nbd in hb:
                num_in += 1
            else:
                num_out += 1
        ratio_nbd[ijk] = float(num_out)/(num_in+num_out)
    return ratio_nbd

def average_nbds(hb, bdry_in, bdry_out, dat, get_nbd):
    factor_bdry = 0.5
    avg_nbd = {}
    for ijk in bdry_in+bdry_out:
        int_in = 0.0
        int_out = 0.0
        num_in = 0.0
        num_out = 0.0
        nbds = get_nbd.get_nbds_26(ijk)
        for nbd in nbds:
            if nbd in bdry_in:
                int_in += dat[nbd] * factor_bdry
                num_in += factor_bdry
            elif nbd in hb:
                int_in += dat[nbd]
                num_in += 1
            elif nbd in bdry_out:
                int_out += dat[nbd] * factor_bdry
                num_out += factor_bdry
            else:
                int_out += dat[nbd]
                num_out += 1
        if ijk in bdry_in:
            int_in += dat[ijk]
            avg_nbd[ijk] = int_in/(num_in+1)
        else:
            int_out += dat[ijk]
            avg_nbd[ijk] = int_out/(num_out+1)
    return avg_nbd
        
def get_partial_volume(hb, dat, get_nbd, verbose):
    partial_volumes = 0.0
    bdry_in = set_bdry(hb, None, get_nbd)
    bdry_out = set_bdry_out(hb, None, get_nbd)
    ratio_nbd = count_nbds(hb, bdry_in, bdry_out, get_nbd)
    avg_nbd = average_nbds(hb, bdry_in, bdry_out, dat, get_nbd)
    pv_dic = {}

    for ijk in hb:
        pv_dic[ijk] = 1.0

    for ijk in bdry_in:
        #partial_volumes -= (1.0 - get_a_partial_volume(ijk, hb, bdry_in, bdry_out, dat, ratio_nbd, get_nbd, verbose) )
        pv_dic[ijk] = get_a_partial_volume(ijk, hb, bdry_in, bdry_out, dat, avg_nbd, get_nbd, verbose, True)
        partial_volumes -= (1.0 - pv_dic[ijk])

    for ijk in bdry_out:
        #partial_volumes += get_a_partial_volume(ijk, hb, bdry_in, bdry_out, dat, ratio_nbd, get_nbd, verbose, True)
        pv_dic[ijk] = get_a_partial_volume(ijk, hb, bdry_in, bdry_out, dat, avg_nbd, get_nbd, verbose, True)
        partial_volumes += pv_dic[ijk]

    return partial_volumes, pv_dic

def apply_lateral_limit(seg_color, hb):
    # X : RL, Y : PA, Z : IS
    # TODO: for general orientation
    pa_index = [ j for (i,j,k) in hb ]

    min_pa = min(pa_index)
    max_pa = max(pa_index)
    hb_cut = []
    for pa in range(min_pa, max_pa+1):
        hb_pa = [ (i,j,k) for (i,j,k) in hb if j == pa ]
        if hb_pa == []:
            continue

        lr_base = [ i for (i,j,k) in hb_pa ]
        if seg_color == 1:
            lr = range( max(lr_base), min(lr_base)-1, -1 )
            lim_lr = -np.inf
        else:
            lr = range( min(lr_base), max(lr_base) + 1 )
            lim_lr = np.inf

        lr_ub = [-np.inf] * len(lr)

        for (i,j,k) in hb_pa:
            index_i = lr.index(i)
            if lr_ub[index_i] < k:
                lr_ub[index_i] = k

        dec = False
        start_i = 0
        while True:
            prev = lr_ub[start_i]
            if np.isfinite(prev):
                break;
            start_i += 1
        
        for i in range(start_i+1, len(lr)):
            if lr_ub[i] < prev:
                dec = True
            elif prev < lr_ub[i] and dec == True:
                lim_lr = lr[i]
                break
            prev = lr_ub[i]

        if np.isfinite(lim_lr):
            if seg_color == 1:
                hb_cut += [ (i,j,k) for (i,j,k) in hb_pa if i > lim_lr ]
            else:
                hb_cut += [ (i,j,k) for (i,j,k) in hb_pa if i < lim_lr ]
        else:
            hb_cut += hb_pa

        #print pa, lr, lr_ub, lim_lr

    return hb_cut

def apply_inferior_limit_lr(opt, gaussians, roi_voxels_lr):
    return  apply_inferior_limit(1, roi_voxels_lr[0], opt.data_t1, opt.data_t2, gaussians[1], gaussians[2], gaussians[4], gaussians[5]), \
            apply_inferior_limit(2, roi_voxels_lr[1], opt.data_t1, opt.data_t2, gaussians[1], gaussians[2], gaussians[4], gaussians[5])

def apply_inferior_limit(seg_color, hb, data_t1, data_t2, avg1, sig1, avg2, sig2):
    # X : RL, Y : PA, Z : IS
    # TODO: for general orientation
    pa_index = [ j for (i,j,k) in hb ]

    min_pa = min(pa_index)
    max_pa = max(pa_index)
    hb_cut = []
    threshold = avg1-2*sig1
    if threshold < 0:
        print 'Wrong avg_T1 and sigma_T1'

    threshold1 = avg1-3*sig1
    if threshold1 < 0:
        threshold1 = threshold / 2

    #print ' CSF T1 threshold: %s (t1mean: %s, t1std: %s)' % (threshold, avg1, sig1)
    for pa in range(min_pa, max_pa+1):
        limits = []
        hb_pa = [ (i,j,k) for (i,j,k) in hb if j == pa ]
        if hb_pa == []:
            continue
        ks_base = [ k for (i,j,k) in hb_pa ]
        ks_base.sort()

        ks = range( max(ks_base), min(ks_base)-1, -1 )
        k_l = [False] * len(ks)
        k_r = [False] * len(ks)
        k_l1 = [False] * len(ks)
        k_r1 = [False] * len(ks)
      
        for (i,j,k) in hb_pa:
            #index_ks = ks.index[k]
            index_ks = ks[0] - k
            if not k_r[index_ks]:
                for lrnbd in [1, 2]:
                    if data_t1[i+lrnbd,j,k] < threshold:    # check T1
                        k_r[index_ks] = True
                        break
            if not k_r1[index_ks]:
                for lrnbd in [1, 2, 3, 4]:
                    if data_t1[i+lrnbd,j,k] < threshold1:    # check T1
                        k_r1[index_ks] = True
                        break
            if not k_l[index_ks]:
                for lrnbd in [-2, -1]:
                    if data_t1[i+lrnbd,j,k] < threshold:    # check T1
                        k_l[index_ks] = True
                        break
            if not k_l1[index_ks]:
                for lrnbd in [-1, -2, -3, -4]:
                    if data_t1[i+lrnbd,j,k] < threshold1:    # check T1
                        k_l1[index_ks] = True
                        break

        for i in range(len(k_r)):
            if k_r1[i] == False:
                k_r[i] = False
            if k_l1[i] == False:
                k_l[i] = False

        #if len([i for i in k_l if i is True]) < len([i for i in k_r if i is True]):
        if seg_color == 1:
            csf_position = np.array(k_r)
        else:
            csf_position = np.array(k_l)

        lim_i = -np.inf
        if csf_position.any():
            jrange = range(csf_position.argmax()+1, len(csf_position))
            count_false = 0
            count_true = 1
            for j in jrange:
                if not csf_position[j]:
                    count_false += 1
                else:
                    count_false = 0
                    count_true += 1
                #if count_false > 2:
                if count_false > 2 and count_true > 1:
                    lim_i = ks[j-2]
                    count_false = 0
                    break
            if count_false > 0:
                lim_i = ks[-count_false]

        if np.isfinite(lim_i):
            limits.append(lim_i)
        #else:
        #    limits.append(ks[-1])

        #print pa, ks, csf_position, lim_i

        if False:
            hb_cut += hb_cut_a
        if True:
            if seg_color == 1:
                csf_position = np.array([ max( [  i for (i,j,k) in hb_pa if k == ks[0] - index_ks] + [-np.inf]) for index_ks in range(len(ks)) ])
            else:
                csf_position = np.array([ max( [ -i for (i,j,k) in hb_pa if k == ks[0] - index_ks] + [-np.inf]) for index_ks in range(len(ks)) ])

            start_dec = False
            #prev = -np.inf
            prev = csf_position.max()
            lim_i_shape = -np.inf
            vertical = 0
            increase = 0
            #print pa, ks, csf_position,

            if (csf_position<np.inf).any() and (csf_position>-np.inf).any():
                jrange = range( (csf_position == csf_position[csf_position<np.inf].max() ).nonzero()[0][0]+1, len(csf_position))
                #print jrange
                for j in jrange:
                    if np.isinf(csf_position[j]):
                        lim_i_shape = ks[j]
                        break
                    if start_dec:
                        if csf_position[j] == prev:
                            if vertical > 2:
                                lim_i_shape = ks[j]
                                break
                            vertical += 1
                        elif csf_position[j] > prev:
                            if csf_position[j] - prev + increase > 1:
                                lim_i_shape = ks[j]
                                break
                            increase += 1
                        else:
                            vertical = 0
                    elif csf_position[j] < prev:
                        start_dec = True
                    prev = csf_position[j]

            if np.isfinite(lim_i_shape):
                limits.append(lim_i_shape-1)
            #else:
            #    limits.append(ks[-1])

        if True:
            k_l = [-np.inf] * len(ks)
            k_r = [np.inf] * len(ks)
            ks = range( max(ks_base), min(ks_base)-1, -1 )
            for (i,j,k) in hb_pa:
                #index_ks = ks.index[k]
                index_ks = ks[0] - k
                for lrnbd in [1, 2, 3, 4]:
                    if data_t1[i+lrnbd,j,k] < threshold:    # check T1
                        if i+lrnbd < k_r[index_ks]:
                            k_r[index_ks] = i+lrnbd
                for lrnbd in [-4, -3, -2, -1]:
                    if data_t1[i+lrnbd,j,k] < threshold:    # check T1
                        if k_l[index_ks] < i+lrnbd:
                            k_l[index_ks] = i+lrnbd
            
            if seg_color == 1:
                csf_position = np.array(k_r)
            else:
                csf_position = np.array([-i for i in k_l])

            start_dec = 0
            prev = csf_position.max()
            lim_i_shape = -np.inf
            vertical = 0
            increase = 0
            #print pa, ks, csf_position,

            if (csf_position<np.inf).any() and (csf_position>-np.inf).any():
                jrange = range( (csf_position == csf_position[csf_position<np.inf].max() ).nonzero()[0][0]+1, len(csf_position))
                #print jrange
                for j in jrange:
                    if np.isinf(csf_position[j]):
                        lim_i_shape = ks[j]
                        break
                    if start_dec > 1:
                        if csf_position[j] == prev:
                            if vertical > 1:
                                lim_i_shape = ks[j]
                                break
                            vertical += 1
                        elif csf_position[j] > prev:
                            if csf_position[j] - prev + increase > 1:
                                lim_i_shape = ks[j]
                                break
                            increase += 1
                        else:
                            vertical = 0

                    if csf_position[j] < prev:
                        start_dec += 1
                    prev = csf_position[j]

            if np.isfinite(lim_i_shape):
                limits.append(lim_i_shape-1)
            #else:
            #    limits.append(ks[-1])

        #print limits
        if True:
            if limits:
#FIXME
                hb_cut += [ (i,j,k) for (i,j,k) in hb_pa if k >= np.mean(limits) ]
                #hb_cut += [ (i,j,k) for (i,j,k) in hb_pa if k >= max(limits) ]
            else:
                hb_cut += hb_pa

    return hb_cut

def center_of_mass(inds):
    center_x = 0.0
    center_y = 0.0
    center_z = 0.0
    for (i,j,k) in inds:
        center_x += i
        center_y += j
        center_z += k
    N = len(inds)
    center_x /= N
    center_y /= N
    center_z /= N
    center_x = int(np.ceil(center_x))
    center_y = int(np.ceil(center_y))
    center_z = int(np.ceil(center_z))
    return center_x, center_y, center_z


def weighted_center_of_mass(inds, dat):
    center_x = 0.0
    center_y = 0.0
    center_z = 0.0
    N = 0.0
    for (i,j,k) in inds:
        w = dat[i,j,k]
        center_x += (i*w)
        center_y += (j*w)
        center_z += (k*w)
        N += w
    center_x /= N
    center_y /= N
    center_z /= N
    center_x = int(np.ceil(center_x))
    center_y = int(np.ceil(center_y))
    center_z = int(np.ceil(center_z))

    #return center_x, center_y, center_z

    def dist(ind):
        return abs(ind[0]-center_x) + abs(ind[1]-center_y) + abs(ind[2]-center_z)

    min_ind = 0
    min_dist = dist(inds[0])
    for i in range(1, len(inds)):
        curr_dist = dist(inds[i])
        if curr_dist < min_dist:
            min_ind = i
            min_dist = curr_dist
        elif curr_dist == min_dist:
            if dat[inds[min_ind]] < dat[inds[i]]:
                min_ind = i
                min_dist = curr_dist

    return inds[min_ind]


def find_max(inds, dat):
    max_ind = 0
    max_val = dat[inds[0]]
    for i in range(1, len(inds)):
        curr_val = dat[inds[i]]
        if max_val < curr_val:
            max_ind = i
            max_val = curr_val

    return inds[max_ind]


class ROI:
    def __init__(self, unit=(0.7,0.7,0.7), dx=3, dy=4, dz=3, template=None, get_nbd=None):

        if template:
            import nibabel as nib
            print 'loading template', template
            nif = nib.load(template)
            dat = nif.get_data()

            seed = (dat == 1).nonzero()
            roi = list(zip(seed[0], seed[1], seed[2]))
            roi = dilate_once(roi, get_nbd)
            center_x, center_y, center_z = center_of_mass(roi)
            self.roi_r = [ (x - center_x, y - center_y, z - center_z) for (x,y,z) in roi ]
            
            seed = (dat == 2).nonzero()
            roi = list(zip(seed[0], seed[1], seed[2]))
            roi = dilate_once(roi, get_nbd)
            center_x, center_y, center_z = center_of_mass(roi)
            self.roi_l = [ (x - center_x, y - center_y, z - center_z) for (x,y,z) in roi ]

        else:
            self.unit = unit
            di = int(dx/unit[0])
            dj = int(dy/unit[1])
            dk = int(dz/unit[2])
            #self.roi_base = [ (i,j,k) for i in range(-di, di+1) for j in range(-dj, dj+1) for k in range(-dk, dk+1)
            #        if dj*dj*dk*dk*i*i + di*di*dk*dk*j*j + di*di*dj*dj*k*k < di*di*dj*dj*dk*dk ]

            # before 4/24/2015
            self.roi_base = [ (i,j,k) for i in range(-di, di+1) for j in range(-dj, dj+1) for k in range(-dk, dk+1)
                    if abs(i) + abs(j) + abs(k) <= max([di,dj,dk]) ]
            roi_r = [ (i+di-1,j,k) for i in range(-di, 1+di/2) for j in range(-dj, dj+1) for k in range(2*dk)
                    if abs(i) + abs(j) + abs(k) <= max([1.5*di,dj,2*dk]) if abs(i) + abs(j) <= max([di, dj]) ]
            roi_l = [ (i-di+1,j,k) for i in range(-di/2, di+1) for j in range(-dj, dj+1) for k in range(2*dk)
                    if abs(i) + abs(j) + abs(k) <= max([1.5*di,dj,2*dk]) if abs(i) + abs(j) <= max([di, dj]) ]

            # template 4/24/2015
            self.roi_base = [ (i,j,k) for i in range(-di, di+1) for j in range(-dj, dj+1) for k in range(-dk, 1)
                    if i*i/float(di*di) + j*j/float(dj*dj) + k*k/float(dk*dk) <= 1.0 ]
            roi_r = [ (i,j,k2+1) for i in range(-di, di+1) for j in range(-dj, dj+1) for k2 in range(2*min([di,dj,dk]))
                    if float((i-k2/2.0)**2)/(di-k2/2.0)**2 + float(j*j)/(dj-k2/2.0)**2 + float(k2/2.0*k2/2.0)/(dk-k2/2.0)**2 <= 1.0 ]
            roi_l = [ (i,j,k2+1) for i in range(-di, di+1) for j in range(-dj, dj+1) for k2 in range(2*min([di,dj,dk]))
                    if float((i+k2/2.0)**2)/(di-k2/2.0)**2 + float(j*j)/(dj-k2/2.0)**2 + float(k2/2.0*k2/2.0)/(dk-k2/2.0)**2 <= 1.0 ]
            

            #roi_l = [ (-i,j,k) for (i,j,k) in roi_r ]
            #roi_r = [ (i+di-1,j,k+dk-1) for i in range(di+1) for j in range(-dj, dj+1) for k in range(-dk, dk+1)
            #        if abs(i) + abs(j) + abs(k) <= max([di,dj,dk])+1 ]
            #roi_l = [ (i-di+1,j,k+dk-1) for i in range(-di,-1) for j in range(-dj, dj+1) for k in range(-dk, dk+1)
            #        if abs(i) + abs(j) + abs(k) <= max([di,dj,dk])+1 ]
            self.roi_r = list(set(self.roi_base + roi_r))
            self.roi_l = list(set(self.roi_base + roi_l))

    def get_roi(self, seed):
        center_x, center_y, center_z = center_of_mass(seed)
        return [ (center_x+i, center_y+j, center_z+k) for (i,j,k) in self.roi_base ]

    def get_roi_lr(self, (center_x, center_y, center_z), lr):
        if lr == 1:
            return [ (center_x+i, center_y+j, center_z+k) for (i,j,k) in self.roi_r ]
        elif lr == 2:
            return [ (center_x+i, center_y+j, center_z+k) for (i,j,k) in self.roi_l ]
        else:
            return False

    def _get_roi_lr(self, (center_x, center_y, center_z), lr):
        if lr == 1:
            return [ (center_x+i + int(1.5/self.unit[0]), center_y+j, center_z+k) for (i,j,k) in self.roi_r ]
        elif lr == 2:
            return [ (center_x+i - int(1.5/self.unit[0]), center_y+j, center_z+k) for (i,j,k) in self.roi_l ]
        else:
            return False

def read_seeds_from_nifti(filename, seed_values=(1,2), verbose=False):
    import nibabel as nib
    nif_s = nib.load(filename)
    dat = nif_s.get_data()
    seed_voxels = []
    for seed_color in seed_values:
        seed = (dat == seed_color).nonzero()
        seed_voxels.append( zip(seed[0], seed[1], seed[2]) )
    if verbose: print '  first seed = (%s, %s) voxels' % (len(seed_voxels[0]), len(seed_voxels[1]))
    return seed_voxels

def segmentation_threshold(opt, roi_voxels_r, roi_voxels_l):
    ''' Return: (thresholded voxels 1, threshold voxels 2), threshold values, gaussian fittings, histograms'''
    roi_voxels = list(set(roi_voxels_r + roi_voxels_l))
    dat1_roi = [ opt.data_t1[ijk] for ijk in roi_voxels ]
    dat2_roi = [ opt.data_t2[ijk] for ijk in roi_voxels ]
    dat3_roi = [ opt.data_my[ijk] for ijk in roi_voxels ]

    ### 1st threshold: exclude CSF

    # Histogram
    bins = 30
    normed = False
    histtype = 'step'
    xx1, y1, a1, avg1, sig1 = histogram_gaussian_fitting(dat1_roi, bins, normed, 'initial T1', fit_min=opt.t1min, verbose=opt.verbose, weight='T1')
    xx2, y2, a2, avg2, sig2 = histogram_gaussian_fitting(dat2_roi, bins, normed, 'initial T2', fit_max=opt.t2max, verbose=opt.verbose, weight='T2')
    xxm, ym, am, avgm, sigm = histogram_gaussian_fitting(dat3_roi, bins, normed, 'initial myelin 1', verbose=opt.verbose)
    # ignore long tail of myelin histogram 
    bins = int(bins/1.5)
    if avgm < 0: rangex = ([avg1/avg2, 4*avg1/avg2])
    else: rangex = (max([avgm-1.5*sigm, 0]), min([max(dat3_roi), avgm+1.5*sigm]))
    xxm, ym, am, avgm, sigm = histogram_gaussian_fitting(dat3_roi, bins, normed, 'initial myelin 2', rangex=rangex, verbose=opt.verbose)

    t1min, t1max, t2min, t2max = set_min_max(avg1, sig1, avg2, sig2, opt.t1min, opt.t1max, opt.t2min, opt.t2max)
    # Threshold
    if opt.verbose: print '  threshold (input or default): t1 in (%s, %s), t2 in (%s, %s)' % (t1min, t1max, t2min, t2max)
    roi_voxels_thr = [ ijk for ijk in roi_voxels if t1min<opt.data_t1[ijk]<t1max and t2min<opt.data_t2[ijk]<t2max ]

    # remove disconnected ROI
    roi_voxels_thr = dilate_intensity(opt.seed_voxels[0]+opt.seed_voxels[1], roi_voxels_thr, opt.data_t1, opt.get_nbd, 0)

    #dat1_roi = [ opt.data_t1[ijk] for ijk in roi_voxels_thr ]
    #dat2_roi = [ opt.data_t2[ijk] for ijk in roi_voxels_thr ]
    dat3_roi = [ opt.data_my[ijk] for ijk in roi_voxels_thr ]

    if opt.verbose:
        print '  Gaussian 1:'
        print '               A      Avg     Sig'
        #                  1234567 1234567 1234567 
        print '    T1      %4.2f %4.2f %4.2f' % (a1, avg1, sig1)
        print '    T2      %4.2f %4.2f %4.2f' % (a2, avg2, sig2)
        print '    myelin  %4.2f %4.2f %4.2f' % (am, avgm, sigm)
        print '  Threshold 1: T1 in (%4.2f, %4.2f), T2 in (%4.2f, %4.2f)' % (t1min, t1max, t2min, t2max)

    ### 2nd threshold: habenula
    # Histogram
    bins = len(roi_voxels_thr)/5
    xx1, y1, a1, avg1, sig1 = histogram_gaussian_fitting(dat1_roi, bins, normed, 'threshold T1', fit_min='auto', verbose=opt.verbose, weight='T1')
    xx2, y2, a2, avg2, sig2 = histogram_gaussian_fitting(dat2_roi, bins, normed, 'threshold T2', fit_max='auto', verbose=opt.verbose, weight='T2')
    xxm, ym, am, avgm, sigm = histogram_gaussian_fitting(dat3_roi, bins, normed, 'threshold myelin 1', verbose=opt.verbose)
    # ignore long tail of myelin histogram 
    bins = int(bins/1.5)
    if avgm < 0:
        rangex = ([avg1/avg2, 4*avg1/avg2])
    else:
        rangex = (max([avgm-1.5*sigm, 0]), min([max(dat3_roi), avgm+1.5*sigm]))
    xxm, ym, am, avgm, sigm = histogram_gaussian_fitting(dat3_roi, bins, normed, 'threshold myelin 2', rangex=rangex, verbose=opt.verbose)

    t1min, t1max, t2min, t2max = set_min_max(avg1, sig1, avg2, sig2, 'avg', opt.t1max, opt.t2min, 'avg')
    # Threshold
    if opt.verbose: print '  threshold (mean of gaussian): t1 in (%s, %s), t2 in (%s, %s)' % (t1min, t1max, t2min, t2max)
    roi_voxels_thr_r = [ ijk for ijk in roi_voxels_r if t1min<opt.data_t1[ijk]<t1max and t2min<opt.data_t2[ijk]<t2max ]
    roi_voxels_thr_l = [ ijk for ijk in roi_voxels_l if t1min<opt.data_t1[ijk]<t1max and t2min<opt.data_t2[ijk]<t2max ]
    # remove disconnected ROI
    roi_voxels_thr_r = dilate_intensity(opt.seed_voxels[0], roi_voxels_thr_r, opt.data_t1, opt.get_nbd, 0)
    roi_voxels_thr_l = dilate_intensity(opt.seed_voxels[1], roi_voxels_thr_l, opt.data_t1, opt.get_nbd, 0)
    # myelin threshold
    min_i = (avg1+opt.sigma_factor*sig1)/(avg2-opt.sigma_factor*sig2)
    roi_voxels_thr_1 = dilate_intensity(opt.seed_voxels[0], roi_voxels_thr_r, opt.data_my, opt.get_nbd, min_i)
    roi_voxels_thr_2 = dilate_intensity(opt.seed_voxels[1], roi_voxels_thr_l, opt.data_my, opt.get_nbd, min_i)

    if opt.verbose:
        print '  Gaussian 2:'
        print '               A      Avg     Sig'
        #                  1234567 1234567 1234567 
        print '    T1      %4.2f %4.2f %4.2f' % (a1, avg1, sig1)
        print '    T2      %4.2f %4.2f %4.2f' % (a2, avg2, sig2)
        print '    myelin  %4.2f %4.2f %4.2f' % (am, avgm, sigm)
        print '  Threshold 2: T1 in (%4.2f, %4.2f), T2 in (%4.2f, %4.2f)' % (t1min, t1max, t2min, t2max)
        print '  Threshold myelin: (%4.2f, )' % min_i
    
    return (roi_voxels_thr_1, roi_voxels_thr_2), (roi_voxels_thr_r, roi_voxels_thr_l), (t1min, t1max, t2min, t2max, min_i), (a1, avg1, sig1, a2, avg2, sig2, am, avgm, sigm), (xx1, y1, xx2, y2, xxm, ym)

def segmentation_region_growth(opt, roi_voxels_thr_lr_myel, roi_voxels_thr_lr_t1t2, up_limit=False):
    pivot_up_limit = 0
    roi_voxels_region_growth_lr = []
    #FIXME
    #data = opt.data_t1
    data = opt.data_my
    for seg_color in [1,2]:
        hb = roi_voxels_thr_lr_myel[seg_color-1][:]
        done = roi_voxels_thr_lr_t1t2[seg_color-1][:]
        for ii in range(opt.region_growth_maxiter):
            hb_extend = dilate_out(hb, done, opt.get_nbd)
            hb_extend = dilate_out(hb_extend, done, opt.get_nbd)

            hb_out = [ijk for ijk in hb_extend if ijk not in hb]
            
            if up_limit and ii == 0:
                pivot_up_limit = calculate_pivot(hb, hb_out, data, 2.0)

            # refine hb (region opt.region_growth_maxiter method)
            stat_prev_mean_i = np.mean([data[ijk] for ijk in hb])
            stat_prev_mean_o = np.mean([data[ijk] for ijk in done if ijk not in hb])

            hb_update, updown = refine_hb(hb, hb_out, data, opt.get_nbd, up_limit=pivot_up_limit, verbose=opt.verbose)

            stat_next_mean_i = np.mean([data[ijk] for ijk in hb_update])
            stat_next_mean_o = np.mean([data[ijk] for ijk in done if ijk not in hb_update])
            if opt.verbose:
                print '  Hb refine: in %s, out %s, mean_in %.2f -> %.2f, mean_out %.2f -> %.2f, mean_diff %.2f -> %.2f' \
                        % (updown[0], updown[1], stat_prev_mean_i, stat_next_mean_i, stat_prev_mean_o, stat_next_mean_o, stat_prev_mean_i-stat_prev_mean_o, stat_next_mean_i-stat_next_mean_o)

            #FIXME
            #hb = hb_update
            if stat_prev_mean_i-stat_prev_mean_o < stat_next_mean_i-stat_next_mean_o:
                hb = hb_update
            else:
                if opt.verbose: print '  Hb refine: Mean differene did not increase. Reject last refine'
                break

            if updown == (0,0) or updown == (1,0) or updown == (0,1):
            #if updown == (0,0):
                break
        
        roi_voxels_region_growth_lr.append(hb)
    return roi_voxels_region_growth_lr

def set_template_roi(opt, roi_voxels_thr_lr, threshold_values):
    roi = ROI(opt.zooms, template=opt.template, get_nbd=opt.get_nbd)
    roi_voxels_thr_t1t2 = []
    roi_voxels_thr_myel = []
    for seg_color in [1,2]:
        #
        # FIXME
        #opt.seed_voxels[seg_color-1] = [weighted_center_of_mass(roi_voxels_thr_lr[seg_color-1], opt.data_my)]

        #roi_template = roi.get_roi_lr(center_of_mass(roi_voxels_thr_lr[seg_color-1]), seg_color)
        roi_template = roi.get_roi_lr(weighted_center_of_mass(roi_voxels_thr_lr[seg_color-1], opt.data_my), seg_color)
        #com = center_of_mass(roi_voxels_thr_lr[seg_color-1])
        #roi_template = dilate_volume([(int(com[0]), int(com[1]), int(com[2]))], opt.get_nbd, opt.min_volume, opt.unit)
        #
        roi_template = [ ijk for ijk in roi_template
                if threshold_values[0]<opt.data_t1[ijk]<threshold_values[1] and threshold_values[2]<opt.data_t2[ijk]<threshold_values[3] ]

        # FIXME
        opt.seed_voxels[seg_color-1] = [find_max(roi_template, opt.data_my)]

        roi_voxels_thr_t1t2.append(dilate_intensity(opt.seed_voxels[seg_color-1], roi_template, opt.data_my, opt.get_nbd, 0))
        roi_voxels_thr_myel.append(dilate_intensity(opt.seed_voxels[seg_color-1], roi_voxels_thr_t1t2[seg_color-1], opt.data_my, opt.get_nbd, threshold_values[4], conn=26))
    return roi_voxels_thr_myel, roi_voxels_thr_t1t2

def segmentation_partial_volume_estimation(opt, roi_voxels_lr, roi_voxels_lr_t1t2):
    roi_voxels_partial_volume_lr = []
    partial_volume_lr = []
    for seg_color in [1, 2]:
        hb_extend = dilate_out(roi_voxels_lr[seg_color-1], roi_voxels_lr_t1t2[seg_color-1], opt.get_nbd)
        hb_extend = dilate_out(hb_extend, roi_voxels_lr_t1t2[seg_color-1], opt.get_nbd)
        hb_out = [ijk for ijk in hb_extend if ijk not in roi_voxels_lr[seg_color-1]]
        partial_volume, pv_dic = get_partial_volume(roi_voxels_lr[seg_color-1], opt.data_my, opt.get_nbd, verbose=opt.verbose)
        roi_voxels_partial_volume_lr.append(pv_dic)
        partial_volume_lr.append(partial_volume * opt.unit)
    return roi_voxels_partial_volume_lr, partial_volume_lr

def make_segmentation(opt, roi_voxels_lr, roi_voxels_partial_volume_lr):
    for seg_color in [1, 2]:
        for ijk in roi_voxels_lr[seg_color-1]:
            opt.out_binary[ijk] = seg_color
        for ijk in roi_voxels_partial_volume_lr[seg_color-1]:
            opt.out_partial[ijk] = roi_voxels_partial_volume_lr[seg_color-1][ijk]
    
def apply_lateral_limit_lr(opt, roi_voxels_lr):
    roi_voxels_constraint_lr = []
    for seg_color in [1, 2]:
        roi_voxels_constraint = apply_lateral_limit(seg_color, roi_voxels_lr[seg_color-1])
        roi_voxels_constraint_lr.append(dilate_intensity(opt.seed_voxels[seg_color-1], roi_voxels_constraint, opt.data_my, opt.get_nbd, 0, conn=26))
    return roi_voxels_constraint_lr

def segmentation(opt, return_volume=True, return_step_volume=False, return_center_of_mass=False, return_hb_index=False, return_hb_neighbor_index=False):
    # using both left and right seed voxels
    if opt.verbose: print '  read initial seed voxels'
    seeds_bi = opt.seed_voxels[0] + opt.seed_voxels[1]
    
    # dilate until > min_volume
    roi_voxels_r = dilate_volume(opt.seed_voxels[0], opt.get_nbd, opt.min_volume, opt.unit)
    roi_voxels_l = dilate_volume(opt.seed_voxels[1], opt.get_nbd, opt.min_volume, opt.unit)
    roi_voxels_bi = list(set(roi_voxels_r + roi_voxels_l))
    if opt.verbose: print '  dilate seeds: %s voxels' % len(roi_voxels_bi)

    # FIXME
    opt.seed_voxels[0] = [weighted_center_of_mass(roi_voxels_r, opt.data_my)]
    opt.seed_voxels[1] = [weighted_center_of_mass(roi_voxels_l, opt.data_my)]

    # segmentation thresholding
    roi_voxels_thr_lr_myel, roi_voxels_thr_lr_t1t2, threshold_values, gaussians, histograms_thr = segmentation_threshold(opt, roi_voxels_r, roi_voxels_l)

    # template ROI
    roi_voxels_thr_lr_myel, roi_voxels_thr_lr_t1t2 = set_template_roi(opt, roi_voxels_thr_lr_myel, threshold_values)
    if opt.verbose: print '  threshold: %s / %s voxels' % (len(roi_voxels_thr_lr_myel[0]), len(roi_voxels_thr_lr_myel[1]))
    if opt.verbose: print '  threshold: %s / %s mm3' % (len(roi_voxels_thr_lr_myel[0])*opt.unit, len(roi_voxels_thr_lr_myel[1])*opt.unit)

    if opt.verbose: print '  threshold_thal: %s / %s voxels' % (len(roi_voxels_thr_lr_t1t2[0])-len(roi_voxels_thr_lr_myel[0]), len(roi_voxels_thr_lr_t1t2[1])-len(roi_voxels_thr_lr_myel[1]))

    # segmentation region growth
    roi_voxels_region_growth_lr = segmentation_region_growth(opt, roi_voxels_thr_lr_myel, roi_voxels_thr_lr_t1t2, up_limit=False)
    if opt.verbose: print '  region growth: %s / %s voxels' % (len(roi_voxels_region_growth_lr[0]), len(roi_voxels_region_growth_lr[1]))
    if opt.verbose: print '  region growth: %s / %s mm3' % (len(roi_voxels_region_growth_lr[0])*opt.unit, len(roi_voxels_region_growth_lr[1])*opt.unit)

    # segmentation inferior limit
    roi_voxels_constraint_lr = apply_inferior_limit_lr(opt, gaussians, roi_voxels_region_growth_lr)
    #roi_voxels_constraint_lr = roi_voxels_region_growth_lr
    if opt.verbose: print '  inferion limit: %s / %s voxels' % (len(roi_voxels_constraint_lr[0]), len(roi_voxels_constraint_lr[1]))
    if opt.verbose: print '  inferion limit: %s / %s mm3' % (len(roi_voxels_constraint_lr[0])*opt.unit, len(roi_voxels_constraint_lr[1])*opt.unit)

    # segmentation lateral limit
    roi_voxels_constraint_lr = apply_lateral_limit_lr(opt, roi_voxels_constraint_lr)
    if opt.verbose: print '  lr limit: %s / %s voxels' % (len(roi_voxels_constraint_lr[0]), len(roi_voxels_constraint_lr[1]))
    if opt.verbose: print '  lr limit: %s / %s mm3' % (len(roi_voxels_constraint_lr[0])*opt.unit, len(roi_voxels_constraint_lr[1])*opt.unit)

    # partial volume estimation
    roi_voxels_partial_volume_lr, partial_volume_lr = segmentation_partial_volume_estimation(opt, roi_voxels_constraint_lr, roi_voxels_thr_lr_t1t2)
    if opt.verbose: print '  partial volume estimation: %s / %s mm3' % (len(roi_voxels_constraint_lr[0])*opt.unit+partial_volume_lr[0], len(roi_voxels_constraint_lr[1])*opt.unit+partial_volume_lr[1])

    # mark segmentation
    make_segmentation(opt, roi_voxels_constraint_lr, roi_voxels_partial_volume_lr)

    return_values = {}
    if return_volume:
        return_values['volumes'] = (len(roi_voxels_constraint_lr[0])*opt.unit+partial_volume_lr[0], len(roi_voxels_constraint_lr[1])*opt.unit+partial_volume_lr[1])

    if return_step_volume:
        return_values['step_volumes'] = (
            (len(roi_voxels_thr_lr_myel[0])*opt.unit, len(roi_voxels_thr_lr_myel[1])*opt.unit),
            (len(roi_voxels_region_growth_lr[0])*opt.unit, len(roi_voxels_region_growth_lr[1])*opt.unit),
            (len(roi_voxels_constraint_lr[0])*opt.unit, len(roi_voxels_constraint_lr[1])*opt.unit),
            (len(roi_voxels_constraint_lr[0])*opt.unit+partial_volume_lr[0], len(roi_voxels_constraint_lr[1])*opt.unit+partial_volume_lr[1])
            )

    if return_center_of_mass:
        return_values['center_of_mass'] = (center_of_mass(roi_voxels_constraint_lr[0]), center_of_mass(roi_voxels_constraint_lr[1]))

    if return_hb_index:
        return_values['hb_index'] = roi_voxels_constraint_lr

    if return_hb_neighbor_index:
        hb_boundary_lr = []
        for seg_color in [1, 2]:
            hb_extend = dilate_out(roi_voxels_constraint_lr[seg_color-1], roi_voxels_thr_lr_t1t2[seg_color-1], opt.get_nbd)
            hb_extend = dilate_out(hb_extend, roi_voxels_thr_lr_t1t2[seg_color-1], opt.get_nbd)
            hb_boundary_lr.append([ijk for ijk in hb_extend if ijk not in roi_voxels_constraint_lr[seg_color-1]])
        return_values['hb_neighbor_index'] = hb_boundary_lr

    return return_values


def segment_hb( t1, t2, mye, seed, out_filename, min_volume=80, t1min=None, t1max=None, t2min=0, t2max=None, sig_factor=0.9, growth=5, template = None, verbose = False, return_volume=True, return_step_volume=False, return_center_of_mass=False, return_hb_index=False, return_hb_neighbor_index=False ):
    '''
        segment_hb( 
            t1,
            t2,
            myelin,
            hb_center,
            out_filename,
            min_volume = 80, # minimum volume of ROI
            t1min = None,
            t1max = None,
            t2min = 10,
            t2max = None,
            sig_factor = 0.9,  # mu + sig_factor*sigma for Hb segmentation
            growth = 5,        # number of iteration for region growth
            verbose = False
            )
        segment Hb.
        input:
            t1, t2, myelin: T1w, T2w, Myelin Nifti1 images
            hb_center: mask Nifti1 file indicating Hb centers (1: right Hb, 2: left Hb)

        output:
            out_filename: output (segmentation) Nifti1 image file name

        parameters:
            min_volume: minimum ROI volume (mm^3) for each Hb. Default 80.
            t1min, t1max, t2min, t2max: first thresold values for T1w, T2w.
            sig_factor: Factor of sigma for myelin threshold for Hb segmentation. Default 1.4
                        The higher factor, the smaller Hb.
            growth: Number of iteration for region growth method. 0 for disable it. Default 2.
            verbose: If True, Display details. Default False.

        return:
            [color1 Hb volume, color2 Hb volume], [color1 Hb partial volume, color2 Hb partial volume], [list(left Hb index), list(right Hb index)]
        '''

    # read nifti1 files
    data_t1, data_t2, data_my, affine, header = read_nifti1_files(t1, t2, mye)
    # read seeds
    seed_voxels = read_seeds_from_nifti(seed, verbose=verbose)
    if len(seed_voxels[0]) == 0 or len(seed_voxels[1]) == 0:
        print 'number of seed voxels is not sufficient'
        return -1

    #FIXME
    #data_t1 = ndimage.generic_filter(data_t1, np.mean, size=3)

    opt = Options(data_t1, data_t2, data_my, seed_voxels, header.get_zooms(),
            min_volume=min_volume, t1min=t1min, t1max=t1max, t2min=t2min, t2max=t2max,
            sigma_factor=sig_factor, region_growth_maxiter=growth, template=template, verbose=verbose)

    return_values = segmentation(opt, return_volume, return_step_volume, return_center_of_mass, return_hb_index, return_hb_neighbor_index)

    # write nifti1 file
    if out_filename[-7:] == '.nii.gz':
        out_filename = out_filename[:-7]
    elif out_filename[-3:] == '.gz':
        out_filename = out_filename[:-3]

    write_nifti1_file(opt.out_binary, out_filename+'.nii.gz', affine, header)
    write_nifti1_file(opt.out_partial, out_filename+'_probability.nii.gz', affine, header)

    return return_values

