import spectrograph as spec
import telescope as tel
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits

def observe_spectrum(instrument, texp, input_wv, input_spec, skyfile="eso_newmoon_radiance.txt"):

    if (os.path.isfile(skyfile)):
        full_path = skyfile
    else:
        full_path = os.environ['COATINGS_PATH']+skyfile
    
    # ESO dark night sky spectrum is in weird units: photons/m2/s/micron/arcsec2
    # interpolate immediately onto the instrument's wavelength grid
    skyspec= np.genfromtxt(full_path,usecols=[0,1],names=['waves_nm','skyflux'])
    sky   = np.interp(instrument.waves,skyspec['waves_nm'],skyspec['skyflux'])

    # Object containing operational parameters of the telescope
    magellan = tel.Telescope()

    # "sky" is in photons/m2/s/micron/arcsec^2, we need to turn this into electrons (e-)
    # which requires multiplying by all factors in the demonimator
    # 
    # (a) magellan.Atel/(100**2) is the telescope area in m2
    # (b) texp is the number of seconds in the exposure
    # (c) instrument.fiber.Afib is the area of one fiber in arcsec^2 (bigger fibers would take in more sky)
    # (d) The tricky one is that "per micron," which is the bandwidth of one pixel. It's not the same as the
    #       dlambda of one pixel, becuase you get light spread from adjacent wavelengths from the wings of the
    #       spectrograph line spread function (LSF).  The spectral resolution R=lambda/dlambda where
    #       dlambda is one "resolution element" and is the effective bandpass of a pixel after convolution with
    #       the LSF. This means that dlambda = lambda / R, and the factoe of 1.0e3 converts from nm to microns.
    #       That's the meaning of the (instrument.waves/1.0e3)/instrument.R line.
    #
    # Finally, to get from sky photons per pixel interval to photoelectrons in the CCD, you need to account
    # for throughput losses in the telescope and instrument.

    skyphotons = sky * \
        magellan.Atel/(100**2) * \
        texp * \
        (instrument.waves/1.0e3)/instrument.R * \
        instrument.fiber.Afib * \
        instrument.throughput * \
        magellan.throughput(instrument.waves) #e-

    # Convert from energy units (ergs) to counted photons
    h = 6.6e-27
    c = 3.0e17   # nm/sec
    input_photons = input_spec / (h * c/input_wv)

    # interpolate the input spectrum onto the LLAMAS wavelength grid
    objspec = np.interp(instrument.waves,input_wv,input_photons)

    # Same calculation as for sky photons, except assume all of the light goes down one fiber
    # So, you don't need to integrate over the fiber area, but do need everything else
    objphotons = objspec * \
        magellan.Atel * \
        texp * \
        (instrument.waves)/instrument.R * \
        instrument.throughput * \
        magellan.throughput(instrument.waves) #e-

    readnoise = instrument.sensor.rn #e-
    dark      = instrument.sensor.dark * texp #e-/s * s

    # This calculates the number of pixels that a fiber subtends. When we simulate the extracted the spectrum, we need
    # to quadrature sum read noise from ALL pixels in the profile, not just one pixel.
    pix_resel = (instrument.fiber.dFib) * (instrument.f_cam / instrument.f_col) / instrument.sensor.pixelsize
    
    skynoise = np.sqrt(skyphotons)
    totnoise = np.sqrt(skyphotons + dark*np.ceil(pix_resel) + readnoise**2*np.ceil(pix_resel))

    return objphotons, totnoise


