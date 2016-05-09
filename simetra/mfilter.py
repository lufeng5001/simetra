import numpy as np

def get_template_average(template, beam):
    """
    Calculate the (beam-weighted) average of the lightcurve template <F>.

    template.shape => (t,)
    beam.shape => (x, y, t)
    output.shape => (x, y)

    where t is the number of time samples, (x, y) is the image dimension. 

    """
    return np.inner(beam, template) / np.sum(beam, axis=-1)

def get_template_dev(template, template_avg):
    """
    Calculate (F - <F>).

    template.shape => (t,)
    template_avg.shape => (x, y)
    output.shape => (x, y, t)

    where t is the number of time samples, (x, y) is the image dimension. 

    """
    # reshape template_avg from (x, y) to (x, y, t)
    return template - np.repeat(template_avg, template.shape[0]).reshape(template_avg.shape+template.shape)

def get_rho(template_dev, data):
    """
    Calculate rho (matched filter statistic).

    template_dev.shape => (x, y, t)
    data.shape => (x, y, t)
    output.shape => (x, y)

    where t is the number of time samples, (x, y) is the image dimension. 

    """
    return np.sum(data*template_dev, axis=-1)

def get_sigma_rho(template_dev, beam):
    """
    Calculate sigma of rho (the width of the rho distribution and also a normalization factor). 

    """
    return np.sqrt(np.sum(beam*template_dev**2, axis=-1))

def run_matched_filter(data, beam, template):
    """
    Return rho and sigma_rho given data, beam, and template. 
    (combine the individual functions in this script)

    Note: data and beam are "weighted" already. 

    """
    template_avg = get_template_average(template, beam)
    template_dev = get_template_dev(template, template_avg)
    rho = get_rho(template_dev, data)
    sigma_rho = get_sigma_rho(template_dev, beam)
    return rho, sigma_rho
