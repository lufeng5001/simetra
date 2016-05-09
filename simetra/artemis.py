import optparse, glob
import simetra as art
import numpy as np
import injection as inj
import logging
from astropy.time import Time

def main(opts):

    # set up logging info
    logging.basicConfig(filename=opts.outfile.split(".")[0]+".log", level=logging.INFO)

    # define image section and assign pixel id
    sec, xid, yid = art.define_image_section(opts.section)
    npzname = opts.section.replace(",","_")+".npz"

    # convert images to pixel lightcurves
    if opts.ioconversion:
        imlist = art.file2list(opts.images)
        imdata, imtime, imnoise = art.timestack(imlist, sec)
        if not opts.beams:
            logging.info("** No primary beam files specified **")
            logging.info("** Analysis will proceed by assuming all primary beam values are 1 **")
            bmdata = np.ones(imdata.shape)
        else:
            bmlist = art.file2list(opts.beams)
            bmdata, notime, nonoise = art.timestack(bmlist, sec, isbeam=True)
        np.savez(npzname, imdata=imdata, imtime=imtime, imnoise=imnoise, bmdata=bmdata)
    else:
        lcarray = np.load(npzname)
        imdata = lcarray["imdata"]
        bmdata = lcarray["bmdata"]
        imnoise = lcarray["imnoise"]
        imtime = Time(lcarray["imtime"])

        # weight the lightcurves
        imdata = imdata * bmdata / imnoise**2  ## ASSUME input images are NOT pbcor images
        bmdata = bmdata**2 / imnoise**2

        # inject transients
        injnpz = None  # default
        if opts.injection: 
            injnpz = opts.injnpz
            logging.info("Inject transients from " + injnpz)
            inj_data, inj_amp, inj_params, inj_t0 = inj.inject_transients(imdata.shape, imtime, injnpz, plot=opts.plot)
            art.writeifits(opts.outfile.split(".")[0]+"_istore.fits", xid, yid, inj_amp, inj_params, inj_t0)
            imdata = imdata + bmdata*inj_data

        # calculate rho
        template_id = map(int, opts.which_template.split(","))  # for single template option
        rho, pnames = art.hunt_transients(imtime, imdata, bmdata, opts.template, opts.single, opts.dtshift, i=template_id[0], j=template_id[1], plot=opts.plot)
        
        # save output
        art.writefits(opts.outfile, xid, yid, rho, pnames, injnpz=injnpz)

if __name__ == "__main__":

    template_choices = ["tophat", "power_law", "fred"]
    parser = optparse.OptionParser()
    parser.set_usage("Usage: python %prog  [options]")
    parser.add_option("-i", "--images", default="", type="str",
                      help="List of input sky images ('*.fits'). [default: %default]")
    parser.add_option("-b", "--beams", default="", type="str",
                      help="List of input primary beam images ('*.fits'). [default: %default]")
    parser.add_option("-s", "--section", default="", type="str",
                      help="Coordinate list to define the image section (x1,x2,y1,y2). [default: %default]")
    parser.add_option("-t", "--template", default="tophat", type="choice", choices=template_choices,
                      help="Template name. Valid options: " + ",".join(template_choices) + " [default: %default]")
    parser.add_option("-o", "--outfile", default="killua.fits", type="str",
                      help="Output filename. Existing files will be overwritten. [default: %default]")
    parser.add_option("-l", "--single", action="store_true", default=False,
                      help="Iterate over one single template. [default: %default]")
    parser.add_option("-c", "--which_template", default="0,0", type="str",
                      help="Single template choice: template_index, start_time_index. [default: %default]")
    parser.add_option("-j", "--injection", action="store_true", default=False,
                      help="Run transient injection. [default: %default]")
    parser.add_option("-f", "--injnpz", default="injparam.npz", type="str",
                      help="Transient injection npz file. [default: %default]")
    parser.add_option("-p", "--plot", action="store_true", default=False,
                      help="Plot search or injection templates. This will slow down the code. [default: %default]")
    parser.add_option("-d", "--dtshift", default=120, type="float",
                      help="Template duration shift in sec. [default: %default (1 snapshot)]")
    parser.add_option("-w", "--ioconversion", action="store_true", default=False,
                      help="Run I/O fits2npz conversion. [default: %default]")
    (opts, args) = parser.parse_args()

    main(opts)
