# Create a pixel mask table that I can combine with artemis output. 
# Do this once for the same set of pixels that I run. 

import optparse
import simetra as art
from astropy.io import fits

def get_pixmask(_pixfile, fitsfile):

    allflags = []
    with open(_pixfile) as pixfile:
        for line in pixfile:
            sec, xid, yid = art.define_image_section(line)  # get pixels that go into artemis
            oneflag, nodate, nonoise = art.loadfits(fitsfile, sec)  # load the masked data
            allflags.extend(list(oneflag.transpose().ravel()))  # put it in the right order
    return allflags

def main(opts, pixfiles):

    flags = []
    for pixfile in pixfiles:
        print "Processing " + pixfile
        flags.extend(get_pixmask(pixfile, opts.infits))
    hdu = fits.new_table([fits.Column(name="flag", array=flags, format="L")])
    hdu.writeto(opts.outfits, clobber=True)

if __name__ == "__main__":

    parser = optparse.OptionParser()
    parser.set_usage("Usage: python %prog [options] *pix.list")
    parser.add_option("-i", "--infits", default="", type="str",
                      help="Input masked pixels fits file. [default: %default]")
    parser.add_option("-o", "--outfits", default="", type="str",
                      help="Output masked pixels fits table. [default: %default]")
    (opts, args) = parser.parse_args()
    main(opts, args)
