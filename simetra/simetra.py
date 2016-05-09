# Collection of functions used in the artemis pipeline.
# simetra = artemis backwards

import sys
from astropy.io import fits
from astropy.time import Time, TimeDelta
import numpy as np
import mfilter
import matplotlib.pyplot as pp
import logging

#################################################################################################################
## INPUT

def define_image_section(secdef):
    """
    Define the image section and assign pixel id. 

    Parameters:
    secdef -- string of the pixel coordinates for the image section (x1, x2, y1, y2)

    Outputs:
    sec -- secdef but converted to a list
    x -- x indices for RA
    y -- y indices for Dec

    """
    sec = map(int, secdef.split(","))  # convert string to list
    logging.info("Processing image section\n\t[ra, dec] = [%i:%i, %i:%i]"%(sec[0],sec[1],sec[2],sec[3]))
    xsize = sec[1] - sec[0]  # image size
    ysize = sec[3] - sec[2]
    x, y = np.indices((xsize, ysize))  # generate index grid
    x = x + sec[0]  # correct for pixel offset, i.e. pixel coordinates not starting at (0,0)
    y = y + sec[2]
    return sec, x, y

def file2list(infile):
    """
    Read input list of fits filenames.

    Parameters:
    infile -- file containing a list of input fits files (in chronological order)

    Outputs:
    flist -- python list of filenames

    """
    flist = []
    with open(infile) as fnames:
        for f in fnames:
            flist.append(f.rstrip())
    return flist

def loadfits(fitsfile, sec):
    """
    Load sectioned image data from an input fits file.

    Parameters:
    fitsfile -- input fits file with 4 naxis (stokes, freq, dec, ra)
    sec -- list of section indices with the format [x1, x2, y1, y2]

    Outputs: 
    fitsdata -- image data for input indices
    obsdate -- time of observation (Time object)
    noise -- estimated noise (Jy)

    """
    ## TODO check stokes/freq axes
    with fits.open(fitsfile) as hdu:
        try:
            fitsdata = hdu[0].section[0, 0, sec[2]:sec[3], sec[0]:sec[1]]  # (stokes, freq, dec, ra)
            #fitsdata = hdu[0].data[0, 0, sec[2]:sec[3], sec[0]:sec[1]]
        except IndexError:
            fitsdata = hdu[0].section[sec[2]:sec[3], sec[0]:sec[1]]  # (dec, ra)
            #fitsdata = hdu[0].data[0, 0, sec[2]:sec[3], sec[0]:sec[1]]
        obsdate = hdu[0].header["DATE-OBS"]
        try:
            noise = hdu[0].header["NOISE"]  # image noise computed elsewhere
        except KeyError:
            noise = np.median(np.absolute(fitsdata.ravel() - np.median(fitsdata.ravel()))) * 1.4826  # mad to std
            #noise = 0.1  # test with wrong noise value
    return fitsdata, obsdate, noise

def timestack(imagelist, sec, isbeam=False):
    """
    Stack sectioned images into a time series. 

    Parameters:
    imagelist -- list of input image file names
    sec -- list of section indices with the format [x1, x2, y1, x2]
    isbeam -- (optional) image list is primary beam

    Outputs:
    farray -- numpy array of the image time series with the shape (xpix, ypix, time)
    tarray -- numpy array of observation times
    narray -- numpy array of noise values (or empty list)

    """
    farray, tarray, narray = [], [], []  # flux array, time array, noise array
    for image in imagelist:
        logging.info("Loading " + image)
        imdata, tobs, noise = loadfits(image, sec)
        farray.append(imdata)
        tarray.append(tobs)
        if not isbeam:
            narray.append(noise)
    return np.array(farray).transpose(), Time(tarray, format="isot", scale="utc"), np.array(narray)

#################################################################################################################

def resample(myarray, shift, include_edge=False):
    """
    Resample an array, e.g. duration or start times, with a certain delta step (shift). 

    Parameters:
    myarray -- input array to resample
    shift -- unit to shift

    Outputs:
    myarray -- resampled

    """
    ishift = int(shift / (myarray[1]-myarray[0]))
    if ishift < 1:  # in case the shift is too small
        ishift = 1
    ix = np.where(np.diff(myarray) > shift)[0] + 1  # indices of gaps
    if include_edge:  ## TODO edge cases for power law where t0 starts before observation (infinite gap left and right)
        ix = np.insert(ix, 0, 0)  # edge indices
        #ix = np.append(ix, myarray.size-1)
    it = np.array([])
    for i in xrange(ix.size):
        try:
            temp = np.arange(ix[i], ix[i+1], ishift)  # index sample
            it = np.append(it, temp)
        except IndexError:
            it = np.append(it, ix[i])  # edge case
    it = it.astype(int)
    return myarray[it]

def loadparams(templatetype, tarray, tshift=600):
    """
    Generate template parameters for a given type of template.

    Parameters:
    templatetype -- input template type
    tarray -- observation times
    tshift -- time step to shift duration (sec)

    Outputs:
    params -- list of tuples sampling the parameter phase space
    pnames -- list of parameter names

    """
    if templatetype == "tophat":
        pnames = ["duration_sec"]  # seconds
        pspace = (tarray[1:] - tarray[0]).sec  # time difference wrt first obs
        pspace = resample(pspace, tshift)
        params = zip(pspace)   ### seconds
    ## TODO add power law phase space
    return params, pnames

def gentemplate(templatetype, param, tarray, t0, plot=False):
    """
    Generate a lightcurve template for a given type of template and a particular set of parameters.

    Parameters:
    templatetype -- input template type
    param -- a tuple of template parameters
    tarray -- observation times for the lightcurve
    t0 -- start time of the transient
    plot -- option to generate plots of the templates (False by default to prevent I/O slowdown)

    Outputs:
    template -- numpy array of template fluxes

    """
    if templatetype == "tophat":
        template = tophat(tarray, t0, param[0])
    if templatetype == "power_law":
        template = power_law(tarray, t0, param[0], param[1])
    if templatetype == "fred":
        template = fred(tarray, t0, param[0], param[1])
    if plot: 
        pp.plot(tarray.mjd, template, marker=".", linestyle="-", color="k")
        pp.xlabel("Time (MJD)")
        pp.ylabel("Flux (Jy)")
        pp.title("Template ("+templatetype+")")
        pp.savefig(templatetype+"_"+str(t0.mjd)+"_"+str(param)+".png")
        pp.clf()
    return template

def hunt_transients(imtime, imdata, bmdata, template_name, single, tshift, i=0, j=0, plot=False):
    """
    Hunt for transients -- the heart of artemis. 
    Two steps:
    1. Generate the template bank.
    2. Maximize the matched filter statistic rho given all template/time combinations.

    Parameters:
    imtime -- numpy array of observation times
    imdata -- image time series
    bmdata -- primary beam time series
    template_name -- template choice
    single -- T/F whether to use a single template (testing purposes)

    Outputs:
    rho -- maximum rho for each pixel
    pnames -- names of the parameters for the template

    """
    logging.info("Generating phase space for template type\n\t" + template_name)
    paramlist, pnames = loadparams(template_name, imtime, tshift)
    if single:  # loop over ONLY one single template and start time
        paramlist = [paramlist[i]]

    logging.info("Calculating rho")
    rhocopy, rhosigcopy = None, None  # holds a copy of rho/significance to compare with the next iteration
    for param in paramlist:  # iterate over template parameters
        logging.info("\tparam = " + str(param))
#        if single:  # loop over ONLY one single template and start time
#            tstart = imtime#[j]  ## all start times, one template
#        else:
        tstart = resample(imtime.gps, param[0]*0.1, include_edge=True)  # sec, frac hardcoded (use 120 for single-snapshot)
        tstart = Time(tstart, format="gps", scale="utc")  # convert tarray back to Time format
        for t in xrange(len(tstart)):  # iterate over start times
            t0 = tstart[t]
            template = gentemplate(template_name, param, imtime, t0, plot=plot)  # generate the template
            rho, sigma_rho = mfilter.run_matched_filter(imdata, bmdata, template)  # compute rho
            rhosig = rho / sigma_rho  # compute rho significance
            logging.info("\t\t"+str(t0.mjd)+"\t"+str(rhosig[0,0]))
            rho = arr2rec(rho, sigma_rho, param, t0)  # convert rho to record array
            if rhocopy is None:  # first iteration
                rhocopy = rho.copy()
                rhosigcopy = rhosig.copy()
            else:  # maximize rho significance
                ix = np.where(rhosig > rhosigcopy)
                rhocopy[ix] = rho[ix]
                rhosigcopy[ix] = rhosig[ix]
    return rhocopy, pnames

def arr2rec(arr, sigma, param, t0):
    """
    Convert the rho array into a record array with the template parameters and the transient time as descriptions.

    Parameters:
    arr -- numpy array (rho)
    sigma -- sigma of rho distribution
    param -- tuple of template parameters
    t0 -- transient time

    Outputs: 
    rho in record array format with comma-separated string for param 'key' (parameters, transient time)

    """
    # convert the parameters into a string and format the string into a workable array
    paramstr = ",".join([str(p) for p in param]) + "," + str(t0.mjd)
    paramarr = np.repeat(np.array(paramstr), arr.size).reshape(arr.shape)
    return np.rec.fromarrays([arr, sigma, paramarr], names="rho, sigma_rho, param")

def writefits(outfile, xpix, ypix, rhoarray, pnames, injnpz):
    """
    Save results to a fits file. 
    Note: pixels are rearranged into 1D because of binary table. 

    Parameters:
    outfile -- output file name
    xpix -- image pixel x index
    ypix -- image pixel y index
    rhoarray -- maximum rho values
    pnames -- parameter names for the template used in the calculation
    injnpz -- injection npz file

    Outputs:
    (none, but outfile is saved to disk)

    """
    # reorganize the template parameters
    parray = []  # parameter array
    for p in np.nditer(rhoarray["param"]):
        # split the parameters from the string bundling them together 
        # this is better for storage and future analyses
        # they were bundled together because rho maximization would be a pain otherwise
        parray.append(map(float, str(p).split(",")))
    parray = np.array(parray)

    # create fits table
    # NOTE xpix and ypix assumes the input was a 2D image, so it might break for 1D arrays
    cols = [fits.Column(name="xpix", array=xpix.ravel(), format="I"),
            fits.Column(name="ypix", array=ypix.ravel(), format="I"),
            fits.Column(name="rho", array=rhoarray["rho"].ravel(), format="E"),
            fits.Column(name="sigma_rho", array=rhoarray["sigma_rho"].ravel(), format="E")]
    for i in xrange(parray.shape[1]-1):
        cols.append(fits.Column(name=pnames[i], array=parray[:,i], format="E"))
    cols.append(fits.Column(name="mjd", array=parray[:,-1], format="D"))  # obstime needs double precision

    # create HDUs
    ## TODO add other header information?
    hdu = fits.new_table(cols)
    if injnpz is not None:
        hdu.header["INJFILE"] = injnpz
    hdu.writeto(outfile, clobber=True)
    return

def writeifits(outfile, xpix, ypix, iamp, iparam, it0):
    """
    Save injection parameters for each pixel. 

    Parameters:
    outfile -- output file name
    xpix -- image pixel x index
    ypix -- image pixel y index
    iamp -- injection amplitude
    iparam -- injection parameters (e.g. duration)
    it0 -- injection start time

    Outputs:
    (none, but outfile is saved to disk)

    """

    # create fits table
    cols = [fits.Column(name="xpix", array=xpix.ravel(), format="I"),
            fits.Column(name="ypix", array=ypix.ravel(), format="I"),
            fits.Column(name="amp", array=iamp.ravel(), format="E"),
            fits.Column(name="param", array=iparam.ravel(), format="E"),
            fits.Column(name="mjd", array=it0.ravel(), format="D")]
    hdu = fits.new_table(cols)
    hdu.writeto(outfile, clobber=True)
    return

#################################################################################################################
## TEMPLATES

def tophat(tarray, t0, dt):
    """
    Generate tophat lightcurve.

    Parameters: 
    tarray -- observation times
    t0 -- start time of the transient
    dt -- duration of the transient

    Outputs:
    lc -- numpy array of template lightcurve

    """
    lc = [0]*len(tarray)  # initialize lightcurve
    t1 = t0 + TimeDelta(dt, format="sec")  # end time of transient
    ix0 = np.where(tarray.jd >= t0.jd)[0][0]  # index of start time
    try:
        ix1 = np.where(tarray.jd >= t1.jd)[0][0]  # index of end time
    except IndexError:  # edge case
        ix1 = len(tarray)
    lc[ix0:ix1] = [1]*(ix1-ix0)  # Jy, same unit as image data
    return np.array(lc)

def power_law(tarray, t0, alpha, beta):
    """
    Generate generic power law (rise and fall) lightcurve. 

    Parameters:
    tarray -- observation times
    t0 -- peak time of the transient
    alpha -- rise index
    beta -- decay index

    Outputs: 
    lc -- numpy array of template lightcurve

    """
    lc = np.zeros(len(tarray))  # initialize lightcurve
    ia = np.where(tarray.jd < t0.jd)[0]  # index before peak time
    ib = np.where(tarray.jd >= t0.jd)[0]  # index after peak time
    lc[ia] = ((tarray[ia]-tarray[0]).jd / (t0-tarray[0]).jd)**alpha  # peak amplitude is 1Jy (same unit as image data)
    lc[ib] = ((tarray[ib]-tarray[0]).jd / (t0-tarray[0]).jd)**beta
    return lc

def fred(tarray, t0, tau1, tau2):
    """
    Generate fast-rise-exponential-decay (FRED) lightcurve, e.g. XRB or M dwarfs. 

    Parameters: 
    tarray -- observation times
    t0 -- start time
    tau1 -- characteristic rise time (in days)
    tau2 -- characteristic decay time (in days)

    Outputs: 
    lc -- numpy array of template lightcurve

    """
    lc = np.zeros(len(tarray))  # initialize lightcurve
    k = np.exp(2*np.sqrt(tau1/tau2))
    ia = np.where(tarray.jd > t0.jd)  # only makes sense for times after start time
    lc[ia] = k * np.exp(-tau1/(tarray[ia].jd-t0.jd) - (tarray[ia].jd-t0.jd)/tau2)  # peak is 1Jy
    return lc
