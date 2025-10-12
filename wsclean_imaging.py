import numpy as np
from casatools import table


def get_freq_from_ms(msname):
    tb = table()
    tb.open(f"{msname}/SPECTRAL_WINDOW")
    chan_freqs = tb.getcol("CHAN_FREQ")             # per-channel frequencies
    tb.close()
    return np.median(chan_freqs.ravel())

def find_smallest_fftw_sz_number(n):
    """
    Find the smallest number that can be decomposed into 2,3,5,7
    
    :param n: input number
    :return: the smallest number that can be decomposed into 2,3,5,7
    """

    max_a = int(np.ceil(np.log(n) / np.log(2)))
    max_b = int(np.ceil(np.log(n) / np.log(3)))
    max_c = int(np.ceil(np.log(n) / np.log(5)))
    max_d = int(np.ceil(np.log(n) / np.log(7)))

    smallest_fftw_sz = float('inf')
    for a in range(max_a + 1):
        for b in range(max_b + 1):
            for c in range(max_c + 1):
                for d in range(max_d + 1):
                    fftw_sz = (2 ** a) * (3 ** b) * (5 ** c) * (7 ** d)
                    if fftw_sz > n and fftw_sz < smallest_fftw_sz:
                        smallest_fftw_sz = int(fftw_sz)
    return smallest_fftw_sz


def make_wsclean_cmd(msfile, imagename, size:int =4096, scale='2arcmin', fast_vis=False, field=None, 
            predict=True, auto_pix_fov = False, telescope_size = 3200, im_fov=200*3600*2/np.pi, pix_scale_factor=1.5,
            **kwargs):  
    """
    Wrapper for imaging using wsclean, use the parameter of wsclean in args. 
    
    To be noted:

    * replace '-' with '_' in the argument name), 
    * The args without value is set to True or False, for example.
    * add False to a key to remove it from the args list.

    ``run_wsclean('OVRO-60MHz.MS', 'IMG-60MHz',  size=2048, niter=1000, mgain=0.9, no_reorder=True, predict=True)``
    
    :param msfile: input CASA measurement set
    :param imagename: output image name
    :param size: (int) image size, default 4096
    :param scale: pixel scale, default 2arcmin
    :param fast_vis: if True, split the measurement set into multiple measurement sets, each containing one field
    :param field: field ID to image, if fast_vis is True
    :param predict: if True, predict the model visibilities from the image
    :param auto_pix_fov: if True, automatically set the pixel scale to match the field of view
    :param telescope_size: size of the telescope in meters, default 3200 (OVRO-LWA)
    :param im_fov: field of view of the image in arcseconds, default 182*3600*2/np.pi (full sky+ 2deg) scaling down by 2/pi
    :param j: number of threads, default 4
    :param mem: fraction of max memory usage, default 2 
    :param weight: weighting scheme, default uniform
    :param no_dirty: don't save dirty image, default True
    :param niter: number of iterations, default 10000
    :param mgain: maximum gain in each cycle, default 0.8
    :param auto_threshold: auto threshold, default 3
    :param auto_mask: auto mask, default 8
    :param pol: polarization, default I
    :param minuv_l: minimum uv distance in lambda, default 10
    :param intervals_out: number of output images, default 1
    :param no_reorder: don't reorder the channels, default True
    """

    
    default_kwargs={
        'j':'16',                    # number of threads
        'mem':'6',                 # fraction of memory usage
        'weight':'uniform',         # weighting scheme
        'no_dirty':'',              # don't save dirty image
        'no_update_model_required':'', # don't update model required
        'no_negative':'',           # no negative gain for CLEAN
        'niter':'10000',            # number of iterations
        'mgain':'0.8',              # maximum gain in each cycle
        'auto_threshold':'3',       # auto threshold
        'auto_mask':'8',            # auto mask
        'pol':'I',                  # polarization
        'minuv_l':'10',             # minimum uv distance in lambda
        'intervals_out':'1',        # number of output images
        'no_reorder':'',            # don't reorder the channels
        'beam_fitting_size':'2',    # beam fitting size
        'horizon_mask':"2deg",      # horizon mask distance (to mask horizon direction RFI)
        'quiet':'',                 # stop printing to stdout, save time
        'save_source_list':'',      # save source list
    }
    
    if fast_vis:
        default_kwargs['weight']='briggs 0.5'

    if auto_pix_fov:
        freq = get_freq_from_ms(msfile)
        scale_num = 1.22*(3e8/freq)/telescope_size * 180/np.pi*3600 / pix_scale_factor
        scale = str(scale_num/60)+'arcmin'
        size = find_smallest_fftw_sz_number(im_fov/scale_num)

    default_kwargs['size']=str(size)+' '+str(size)
    default_kwargs['scale']=scale
    # remove the key if val is False from kwargs
    for key, value in kwargs.items():
        if value is False:
            default_kwargs.pop(key, None)
        elif value is True:
            # Add the key with an empty string as value if True
            default_kwargs[key] = ''
        else:
            default_kwargs[key] = str(value)

    if fast_vis==True:
        if field is None:
            default_kwargs['intervals_out']='1'
            default_kwargs['field']='all'
        else:
            default_kwargs["intervals_out"] =str(len(field.split(',')))
            default_kwargs['field']='all' # magic, has to be 'all', otherwise only 1st time slot has image
    else:
        default_kwargs['intervals_out']='1'
        default_kwargs['field']='all'
    if default_kwargs['intervals_out']!='1' and predict:
        raise RuntimeError("Prediction cannot be done with multiple images.")
    
    cmd_clean = "wsclean "
    # Add additional arguments from default_params
    for key, value in default_kwargs.items():
        # Convert Python-style arguments to command line format
        cli_arg = key.replace('_', '-')
        cmd_clean += f" -{cli_arg} {value} " if value != '' else f" -{cli_arg} "

    if ('I' in default_kwargs['pol']) and ('join_polarizations' not in default_kwargs.keys()) and \
            ('Q' in default_kwargs['pol'] or 'U' in default_kwargs['pol'] \
            or 'V' in default_kwargs['pol']):                                                      
          cmd_clean+= " -join-polarizations "                                                        
    elif (default_kwargs['pol']=='I' or default_kwargs['pol']=='XX' or default_kwargs['pol']=='YY'
              or default_kwargs['pol']=='XX,YY') and ('no_negative' not in default_kwargs.keys()):  
          cmd_clean+= " -no-negative "                                                              

    cmd_clean += " -name " + imagename + " " 
    
    return cmd_clean