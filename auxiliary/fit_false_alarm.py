# This script determines the false alarm probability of rho. 

import optparse
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as pp
import scipy.special
from scipy.optimize import fmin

def negloglikelihood(params, func, x, y):
    log_mu = func(params, x)  # return log of the mean function
    return -np.sum(y*log_mu - np.exp(log_mu) - scipy.special.gammaln(y+1))  # Poisson regression

def log_exponential(params, x):
    return np.log(params[0]) - (x/params[1])#**params[2]  # exponential fit
    #return np.log(norm) - np.log(rhostar) - scipy.special.gammaln(alpha+1) + alpha*np.log(rhoratio) - rhoratio

def norm_cdf(y, n):
    return np.exp(np.log(y)-np.log(y)[-1]) * n  # normalize background cdf to the search region

def get_rho_threshold(pfa, params):  # for the log_exponential function
    return params[1] * (np.log(params[0]) - np.log(pfa))#**(1/params[2])  # rho that corresponds to prob of false alarm pfa

def get_beta_function(a, b):
    return np.exp(scipy.special.gammaln(a)+scipy.special.gammaln(b)-scipy.special.gammaln(a+b))

def get_error(cdf):  # eq. 34 of efficiency error
    b1 = get_beta_function(cdf+1, cdf[-1]-cdf+1)
    b2 = get_beta_function(cdf+2, cdf[-1]-cdf+1)
    b3 = get_beta_function(cdf+3, cdf[-1]-cdf+1)
    return np.sqrt(b3/b1 - (b2/b1)**2) * cdf[-1]

def rho2flux(rho, sigma):
    return rho / sigma * 1000  # Jy to mJy (rho already divided by a sigma)

def main(opts):

    # get rho
    with fits.open(opts.infile) as infile:
        data = infile[1].data
    ix = np.where(data["flag"] == False)  # mask sources
    x, y = data["xpix"][ix], data["ypix"][ix]  # x and y pixel coordinates
    rhosigma = data["sigma_rho"][ix]
    if opts.significance:
        rho = data["rho"][ix]/rhosigma  # normalized
    else:
        rho = data["rho"][ix]

    # playground data
    if opts.pgfile is not None:
        with fits.open(opts.pgfile) as pgfile:
            pgdata = pgfile[1].data
        iy = np.where(pgdata["flag"] == False)
        pgrho = pgdata["rho"][iy] / pgdata["sigma_rho"][iy]
        pgrho_x = np.sort(pgrho)[::-1]
        if opts.drop is not None:  # drop the rhos above threshold
            idrop = np.where(pgrho_x < opts.drop)
            pgrho_x = pgrho_x[idrop]
        pgcdf_y = (np.arange(pgrho_x.size)+1)/opts.nbeams

    # calculate cumulative distribution of rho
    rho_x = np.sort(rho)[::-1]  # CDF from the right
    cdf_y = np.arange(rho.size)+1  # empirical CDF "without" binning (it's binning by the data)
    cdf_y = norm_cdf(cdf_y, opts.npix/opts.nbeams)  # normalize by the search region and synthesized beam size
    cdf_std = np.sqrt(cdf_y)
    negreg = cdf_y - cdf_std  # replace negative values with small, almost zero values
    ix = np.where(negreg <= 0)
    negreg[ix] = 1e-6
    if opts.pgfile is not None:
#        pgstd = np.sqrt(pgcdf_y)
#        pgnegreg = pgcdf_y - pgstd
#        pgix = np.where(pgnegreg <= 0)
#        pgnegreg[pgix] = 1e-6
        print pgcdf_y[-1], cdf_y[-1]
    if opts.rho is not None: 
        print pgrho_x[:opts.rho]

    # fit exponential function
    params_guess = [rho_x.size, 4]#, 1]  # initial guess for the fit parameters (A, rho*, gamma)
    fitparams = fmin(negloglikelihood, params_guess, args=(log_exponential, rho_x[:opts.tail], cdf_y[:opts.tail]), maxfun=1e5, maxiter=1e5)

    # get threshold
    rho_threshold = get_rho_threshold(opts.pfa, fitparams)
    flux_sensitivity = rho2flux(rho_threshold, rhosigma)  # calculate flux sensitivity from median(flux)... rhosigma distribution is somewhat uniform...

    # plot false alarm probability
    pp.figure(figsize=(11,8.25))
    pp.rc("text", usetex=True)
    pp.plot(rho_x, cdf_y, "k-", label="playground", linewidth=2)
    #pp.fill_between(rho_x[:opts.tail], negreg[:opts.tail], cdf_y[:opts.tail] + cdf_std[:opts.tail], facecolor="y", alpha=0.5)
    pp.axvline(rho_x[opts.tail], linestyle=":", color="k")  # line to mark the part of the data I used for fitting
    if opts.pgfile is not None:
       pp.plot(pgrho_x[16:], pgcdf_y[16:], color="r", linestyle="-", linewidth=2)
       pp.plot(pgrho_x[:16], pgcdf_y[:16], color="r", marker=".", label="search", linestyle="none", linewidth=1)
       #pp.fill_between(pgrho_x[:opts.tail], pgnegreg[:opts.tail], pgcdf_y[:opts.tail] + pgstd[:opts.tail], facecolor="y", alpha=0.5)
    if opts.significance:
        pp.xlabel(r"$\tilde{\rho}$", fontsize=16)
    else:
        pp.xlabel(r"$\rho$")
    pp.ylabel(r"$N(\ge \tilde{\rho})$", fontsize=16)
    if opts.plotlog:
        if not opts.significance:
            pp.xscale("log")  # log rho only if not-significance
        pp.yscale("log")  # symlog to allow for negative values
    fit_x = np.linspace(0,13,50)
    fit_y = np.exp(log_exponential(fitparams, fit_x))
    fitstd = np.sqrt(fit_y)
    fitnegreg = fit_y - fitstd
    fitix = np.where(fitnegreg <= 0)
    fitnegreg[fitix] = 1e-6
    pp.plot(fit_x, fit_y, "b--", label="fit", linewidth=2)
    pp.fill_between(fit_x[:opts.tail], fitnegreg[:opts.tail], fit_y[:opts.tail] + fitstd[:opts.tail], facecolor="y", alpha=0.3)
    pp.title(r"$(A, \rho^{\ast}, \gamma) = $" + str(fitparams))
    pp.suptitle(r"For $N(>\rho) = Ae^{-(\rho/\rho^{\ast})^{\gamma}}$")
    pp.annotate(r"$\rho_{th} = $ "+str(round(rho_threshold, 2)), xy=(0.7,0.8), xycoords="figure fraction", fontsize=20)
    pp.annotate(r"median $S_{th} = $ "+str(round(np.median(flux_sensitivity),2)) + " mJy", xy=(0.5,0.7), xycoords="figure fraction", fontsize=20)
    #pp.legend(loc=1)
    pp.ylim([1e-3,1e4])
    pp.tick_params(axis="both", which="both", labelsize=16)
    pp.savefig(opts.outim)
    pp.clf()

    print "Flux threshold = " + str(round(np.median(flux_sensitivity),2)) + " mJy +/- " + str(round(np.median(np.absolute(flux_sensitivity - np.median(flux_sensitivity))) * 1.4826, 2))  # error calculation may be wrong


if __name__ == "__main__":

    parser = optparse.OptionParser()
    parser.set_usage("Usage: python %prog [options]")
    parser.add_option("-i", "--infile", default="", type="str",
                      help="Input artemis file. [default: %default]")
    parser.add_option("-s", "--significance", action="store_true", default=False,
                      help="Use rho significance instead. [default: %default]")
    parser.add_option("-l", "--plotlog", action="store_true", default=False,
                      help="Plot log axis. [default: %default]")
    parser.add_option("-o", "--outim", default="rho_cdf.png", type="str",
                      help="Output filename for CDF plot. [default: %default]")
    parser.add_option("-p", "--pfa", default=1e-3, type="float",
                      help="Probability of false alarm for rho. [default: %default]")
    parser.add_option("-b", "--nbeams", default=16, type="float",
                      help="Number of pixels in each synthesized beam. [default: %default]")
    parser.add_option("-n", "--npix", default=1e6, type="float",
                      help="Number of pixels in the search region. [default: %default]")
    parser.add_option("-t", "--tail", default=10, type="int",
                      help="Number of tail elements to include for the fit. [default: %default]")
    parser.add_option("-a", "--pgfile", default=None, type="str",
                      help="Playground artemis file. [default: %default]")
    parser.add_option("-r", "--rho", default=None, type="int",
                      help="If this is a search run, print the top r rhos. [default: %default]")
    parser.add_option("-d", "--drop", default=None, type="float",
                      help="Drop the rhos above playground/search threshold (input threshold). [default: %default]")
    (opts, args) = parser.parse_args()
    main(opts)
