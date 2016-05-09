# Join the many artemis output tables into one and plot the pixels sampled in the image. 

import optparse
from astropy.io import fits
import matplotlib.pyplot as pp

def jointb(opts, args):
    flist = open(args[0])  # input file list
    f1 = flist.readline().rstrip()  # first file
    print "Appending " + f1
    hdu1 = fits.open(f1)
    tb1 = hdu1[1]
    for f2 in flist:  # subsequent files
        f2 = f2.rstrip()
        print "Appending " + f2
        hdu2 = fits.open(f2)
        tb2 = hdu2[1]
        nrows1 = tb1.header["NAXIS2"]
        nrows2 = tb2.header["NAXIS2"]
        nrows = nrows1 + nrows2
        tb = fits.new_table(tb1.columns, nrows=nrows)
        for name in tb1.columns.names:  # append tables
            tb.data.field(name)[nrows1:] = tb2.data.field(name)
        hdu1.close()
        hdu2.close()
        tb1 = tb
    tb.writeto(opts.outfile, clobber=True)
    flist.close()
    return

def plotpix(opts):
    outim = opts.outfile.split(".")[0]
    with fits.open(opts.outfile) as tb:
        data = tb[1].data
        pp.plot(data['xpix'], data['ypix'], 'k.')  # plot the pixels by index
        pp.xlabel("RA pixel")
        pp.ylabel("Dec pixel")
        pp.savefig(outim+"_zoom.png")
        pp.xlim((0,opts.xdim))
        pp.ylim((0,opts.ydim))
        pp.savefig(outim+".png")
        pp.clf()
    return

def main(opts, args):
    jointb(opts, args)
    plotpix(opts)
    return

if __name__ == "__main__":

    parser = optparse.OptionParser()
    parser.set_usage("Usage: python %prog [options] image.list")
    parser.add_option("-o", "--outfile", default="madao.fits", type="str",
                      help="Output file with input fits tables appended into one table. [default: %default]")
    parser.add_option("-x", "--xdim", default=4096, type="int",
                      help="Image x-dimension in pixels. [default: %default]")
    parser.add_option("-y", "--ydim", default=4096, type="int",
                      help="Image y-dimension in pixels. [default: %default]")
    (opts, args) = parser.parse_args()
    main(opts, args)
