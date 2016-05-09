# Add a mask column to artemis data file. 

import optparse
from astropy.io import fits

def main(opts):
    with fits.open(opts.infile) as artable, fits.open(opts.maskfile) as masktable:
        newcols = artable[1].columns + masktable[1].columns
        hdu = fits.new_table(newcols)
        hdu.writeto(opts.outfile, clobber=True)

if __name__ == "__main__":

    parser = optparse.OptionParser()
    parser.set_usage("Usage: python %prog [options]")
    parser.add_option("-i", "--infile", default="", type="str",
                      help="Input artemis file. [default: %default]")
    parser.add_option("-m", "--maskfile", default="", type="str",
                      help="Input mask file. [default: %default]")
    parser.add_option("-o", "--outfile", default="", type="str",
                      help="Output artemis file. [default: %default]")
    (opts, args) = parser.parse_args()
    main(opts)
