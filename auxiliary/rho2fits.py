# This script converts rho table into fits image. 

import optparse
from astropy.io import fits
import numpy as np

def main(opts):
    
    with fits.open(opts.infile) as artemis:
        data = artemis[1].data
    rhosig = data["rho"]/data["sigma_rho"]
    if opts.flag:
        ix = np.where(data["flag"] == True)
        rhosig[ix] = 0
    
    with fits.open(opts.outfits) as im:
        imdata = im[0].data
    imdata[:] = 0
    imdata[0,0,data["ypix"],data["xpix"]] = rhosig
    with fits.open(opts.outfits, mode="update") as im:
        im[0].data = imdata
        im.flush()

if __name__ == "__main__":

    parser = optparse.OptionParser()
    parser.set_usage("Usage: python %prog [options]")
    parser.add_option("-i", "--infile", default="", type="str",
                      help="Input artemis file. [default: %default]")
    parser.add_option("-o", "--outfits", default="", type="str",
                      help="Output fits image. [default: %default]")
    parser.add_option("-f", "--flag", action="store_true", default=False,
                      help="Set flagged pixels to 0. [default: %default]")
    (opts, args) = parser.parse_args()
    main(opts)
