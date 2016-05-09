# main: generate the injection file with template amplitude, start time, and parameters specific to each template shape
# inject_transients: generate the injected transient lightcurves for artemis.py 

import optparse
from astropy.io import fits
from astropy.time import Time
from astropy.time import TimeDelta
import numpy as np
import simetra as art

def inject_transients(data_shape, tarray, injnpz, plot=False):
    """
    Create injected transient lightcurves. 

    Input: 
    data_shape - shape of the imdata array
    tarray - time array for transient lightcurve
    injnpz - injection npz file

    """
    inj = np.zeros(data_shape)  # initialize injection data array
    inj_params = np.load(injnpz)  # get injection file parameters
    xpix = data_shape[0] / inj_params["amp"].size  # split the pixels into n chunks for n templates (in the x direction)
    amp = np.zeros((data_shape[0], data_shape[1]))  # initialize injection parameter storage arrays
    params = np.zeros((data_shape[0], data_shape[1]))
    t0 = np.zeros((data_shape[0], data_shape[1]))
    for i in xrange(inj_params["amp"].size):
        inj_template = art.gentemplate(inj_params["templatetype"][0], inj_params["params"][i], tarray, inj_params["t0"][i], plot=plot)
        inj[xpix*i:xpix*(i+1),:,:] = inj_template * inj_params["amp"][i]
        amp[xpix*i:xpix*(i+1),:] = inj_params["amp"][i]  # store injection parameters for these pixels
        params[xpix*i:xpix*(i+1),:] = inj_params["params"][i]#[1]  # TEMP too lazy to deal with 2D param array
        t0[xpix*i:xpix*(i+1),:] = inj_params["t0"][i].mjd
    return inj, amp, params, t0

#########################################################################################
## functions below generate the injection file

def loadtime(im):
    """
    Load the header time of the image in Time format. 

    """
    with fits.open(im) as hdu:
        obsdate = hdu[0].header["DATE-OBS"]
    return Time(obsdate, format="isot", scale="utc")

def draw_amp(maxamp, ntemplates, distribution, inarray=[1]):
    """
    Draw random amplitudes from a particular distribution. 

    Inputs: 
    maxamp - maximum amplitude in Jy
    ntemplates - number of draws (for the number of templates)

    """
    if distribution == "special":
        amp_jy = np.array(inarray)  # will be multiplied by maxamp in Jy
    if distribution == "uniform":
        amp_jy = np.random.random_sample(ntemplates)  # 0 to 1 Jy
    return amp_jy * maxamp  # uniform from 0 to maxamp Jy

def draw_t0(tmin, tmax, ntemplates, distribution, inarray=[0]):
    """
    Draw random start times for the transient from a particular distribution.

    Inputs: 
    tmin - time of the first image
    tmax - time of the last image
    ntemplates - number of draws (for the number of templates)

    """
    if distribution == "special":
        trel = np.array(inarray)  # will give tstart at first image
    if distribution == "uniform":
        trel = np.random.random_sample(ntemplates)  # 0 to 1
        #npz = np.load("/nfs/mwa-03/d1/lufeng/pixy/fullset/clean/180MHz_XX_clean_v2/2994_3080_1532_1618.npz")
        #all_t0 = npz["imtime"]
        #ix = np.random.randint(all_t0.size, size=ntemplates)
        #tstart = all_t0[ix]
    tstart = (tmax - tmin) * trel + tmin  # convert to absolute time
    return tstart

def draw_params(tmin, tmax, templatetype, ntemplates, distribution, inarray=[[0.5]]):
    """
    Draw random parameters for a particular template from particular distributions. 

    Inputs:
    tmin - start time of the transient 
           (not time of the first image, because you can't tell if the transient duration runs beyond the last image)
    tmax - time of the last image
    templatetype - the template name
    ntemplates - number of draws (for the number of templates)

    """
    if templatetype == "tophat":
        if distribution == "special":
            tdur = np.array(inarray[0])  # in seconds
        if distribution == "uniform":
            tdur = np.random.random_sample(ntemplates) * (tmax - tmin).sec  # scale to max duration in sec
            #tdur = np.random.random_sample(ntemplates) * 600. + 120.
            #tdur = np.random.random_sample(ntemplates) * 3600. + 3600.
            #tdur = np.random.random_sample(ntemplates) * 7689600.+86400.  # scale to max duration in sec
        params = zip(tdur)
    if templatetype == "power_law":
        if distribution == "special":
            alpha = np.array(inarray[0])  # default: 0.5
            beta = np.array(inarray[1])  # default: -1
        if distribution == "uniform":
            alpha = np.random.random_sample(ntemplates)
            beta = np.random.random_sample(ntemplates)*-1 - 1  # uniform between (-1,-2)
        params = zip(alpha, beta)
    if templatetype == "fred":
        if distribution == "special":
            tau1 = np.array(inarray[0])  ## TODO not sure this works but ok
            tau2 = np.array(inarray[1])
        if distribution == "uniform":
            ## TODO temporarily hardcoded: tau1 ~ 1-2 days, tau2 ~ 30-40 days
            tau1 = np.random.random_sample(ntemplates) * 1. + 1.  # days    (multiply by range then offset by minimum)
            tau2 = np.random.random_sample(ntemplates) * 10. + 30.  # days
        params = zip(tau1, tau2)
    return params

def main(opts):

    # get time range
    images = opts.infile.split(",")
    tstart = loadtime(images[0])
    tend = loadtime(images[1])

    # adjust start time
    if opts.template == "fred":
        tstart = tstart - TimeDelta(864000, format="sec")  # temp: -10 days

    # convert string to lists
    amp_array = map(float, opts.amp.split(","))
    t0_array = map(float, opts.t0.split(","))
    params_array = []
    for parray in opts.param.split(";"):
        params_array.append(map(float, parray.split(",")))
    if len(amp_array) != len(t0_array) or len(amp_array) != len(params_array[0]):
        print "="*50
        print "Warning: You don't have the same number of elements for each template!!!!" 
        print "Ritz is a lazy cheesecake so you need to make sure they're the same otherwise artemis will break (probably)." 
        print "Just copy-pasta the same input like 0,0,0,0 because if you're using special, you're probably not inputting many anyway." 

    # draw template parameters
    amp = draw_amp(opts.maxamp, opts.ntemplates, opts.amp_distribution, inarray=amp_array)  # Jy
    t0 = draw_t0(tstart, tend, opts.ntemplates, opts.t0_distribution, inarray=t0_array)  # whole time range available
    params = draw_params(tstart, tend, opts.template, opts.ntemplates, opts.param_distribution, inarray=params_array)  # subset time range starting from t0

    # save transient parameters
    np.savez(opts.outnpz, amp=amp, t0=t0, params=params, templatetype=[opts.template], distribution=[opts.amp_distribution, opts.t0_distribution, opts.param_distribution])

if __name__ == "__main__":

    template_choices = ["tophat", "power_law", "fred"]
    distribution_choices = ["special", "uniform"]
    parser = optparse.OptionParser()
    parser.set_usage("Usage: python %prog [options]")
    parser.add_option("-i", "--infile", default="im0,im1", type="str",
                      help="Input files: first image and last image. [default: %default]")
    parser.add_option("-o", "--outnpz", default="injparam.npz", type="str",
                      help="Output npz file. [default: %default]")
    parser.add_option("-t", "--template", default="tophat", type="choice", choices=template_choices,
                      help="Template name. Valid options: " + ",".join(template_choices) + " [default: %default]")
    parser.add_option("-n", "--ntemplates", default=1, type="int",
                      help="Number of templates. [default: %default]")
    parser.add_option("-m", "--maxamp", default=0.1, type="float",
                      help="Max amplitude (template is in Jy). [default: %default]")
    parser.add_option("-a", "--amp_distribution", default="uniform", type="choice", choices=distribution_choices,
                      help="Distribution for random amplitude drawing. Valid options: " + ",".join(distribution_choices) + " [default: %default]")
    parser.add_option("-b", "--amp", default="1", type="str",
                      help="Input list of amplitude (multiples of max amplitude). [default: %default]")
    parser.add_option("-d", "--t0_distribution", default="uniform", type="choice", choices=distribution_choices,
                      help="Distribution for random t0 drawing. Valid options: " + ",".join(distribution_choices) + " [default: %default]")
    parser.add_option("-e", "--t0", default="0", type="str",
                      help="Input list of t0 (relative to first image, between 0 and 1). [default: %default]")
    parser.add_option("-p", "--param_distribution", default="uniform", type="choice", choices=distribution_choices,
                      help="Distribution for random parameter drawing. Valid options: " + ",".join(distribution_choices) + " [default: %default]")
    parser.add_option("-q", "--param", default="0.5", type="str",
                      help="Input list of parameters; separate different parameters by semi-colon. [default: %default]")
    (opts, args) = parser.parse_args()
    main(opts)
