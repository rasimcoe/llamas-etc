import spectrograph as spec
import telescope as tel
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits


def enclosed_energy(fwhm_arcsec, r_arcsec, psf='moffat', beta=2.5):
    """Fraction of a point source's light within a circular aperture of radius ``r_arcsec`` for a
    seeing PSF of the given FWHM. ``psf='gaussian'`` or ``'moffat'`` (``beta`` = Moffat power index).
    Returns 1.0 for zero/negative FWHM (perfect seeing -> all light enclosed)."""
    if fwhm_arcsec is None or fwhm_arcsec <= 0:
        return 1.0
    if psf == 'gaussian':
        sigma = fwhm_arcsec / 2.35482
        return 1.0 - np.exp(-(r_arcsec ** 2) / (2.0 * sigma ** 2))
    alpha = fwhm_arcsec / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))       # Moffat
    return 1.0 - (1.0 + (r_arcsec / alpha) ** 2) ** (1.0 - beta)


def aperture_radius(nfib, Afib):
    """Circular-equivalent radius [arcsec] enclosing ``nfib`` fibres of solid angle ``Afib`` [arcsec^2]."""
    return np.sqrt(nfib * Afib / np.pi)


def observe_spectrum(instrument, texp, input_wv, input_spec, airmass=None,
                     source='point', seeing=None, aperture='optimal', psf='moffat', beta=2.5,
                     nbin=1, skyfile="eso_newmoon_radiance.txt", extfile="lco_extinction.txt"):
    """Simulate a LLAMAS observation; return (counts, noise) in e- on the channel's wavelength grid.

    source='point' (default): ``input_spec`` is flux density f_lambda [erg/cm2/s/A]. The fraction of the
        source captured by the fibre aperture is set by the seeing PSF (enclosed energy), and read/dark/sky
        noise is summed over the fibres in the aperture.
          seeing   : FWHM [arcsec]; defaults to 0.8" with a warning.
          aperture : 'single' (1 fibre) | 'optimal' (maximise band SNR) | int (fixed N fibres).
          psf,beta : seeing profile, 'moffat' (default, power ``beta``) or 'gaussian'.
        (To reproduce the pre-v1.1 "all light in one fibre" behaviour, pass seeing=0.)

    source='extended': ``input_spec`` is surface brightness [erg/cm2/s/A/arcsec2]; per-fibre signal =
        SB * fibre area, with no aperture loss. ``nbin`` co-adds nbin fibres (SNR improves as sqrt(nbin)).
    """
    # --- airmass / atmospheric extinction (v1.0) ---
    if airmass is None:
        airmass = 1.0
        print("WARNING: no airmass given -- assuming airmass = 1.0 (zenith, best case). "
              "Pass airmass=sec(z) for your observation; higher airmass gives lower throughput.")
    else:
        print("   applying atmospheric extinction at airmass = {:.3f}".format(airmass))

    # --- sky spectrum (ESO dark-night radiance, photons/m2/s/micron/arcsec2) ---
    if (os.path.isfile(skyfile)):
        full_path = skyfile
    else:
        try:
            coat_path = os.environ['COATINGS_PATH']
        except:
            os.environ['COATINGS_PATH'] = './COATINGS'
        full_path = os.environ['COATINGS_PATH'] + skyfile
    skyspec = np.genfromtxt(full_path, usecols=[0, 1], names=['waves_nm', 'skyflux'])
    sky = np.interp(instrument.waves, skyspec['waves_nm'], skyspec['skyflux'])

    magellan = tel.Telescope()

    # Telescope mirror factor: gated to 1 in measured mode (already in the curve); applied in theoretical.
    if getattr(instrument, 'throughput_mode', 'theoretical') == 'measured':
        tel_throughput = np.ones_like(instrument.waves)
    else:
        tel_throughput = magellan.throughput(instrument.waves)

    # Atmospheric transmission for the source: 10^(-0.4 k(lambda) airmass).
    ext_path = extfile if os.path.isfile(extfile) else os.environ.get('COATINGS_PATH', './COATINGS/') + extfile
    try:
        ke = np.genfromtxt(ext_path, usecols=[0, 1], names=['wave_nm', 'k'])
        kwave = np.interp(instrument.waves, ke['wave_nm'], ke['k'], left=ke['k'][0], right=ke['k'][-1])
        atm = 10.0 ** (-0.4 * kwave * airmass)
    except Exception as e:
        print("   extinction file unavailable (" + str(e) + "); NOT applying atmospheric extinction")
        atm = np.ones_like(instrument.waves)

    # Sky electrons in ONE fibre (per resolution element). See v1.0 notes on the (waves/1e3)/R bandwidth.
    sky_perfib = sky * magellan.Atel / (100 ** 2) * texp * (instrument.waves / 1.0e3) / instrument.R * \
        instrument.fiber.Afib * instrument.throughput * tel_throughput

    # Source: convert f_lambda (or surface brightness) to photons, interpolate onto the grid, and carry
    # through area * texp * bandwidth * throughput * telescope * extinction.
    h = 6.6e-27
    c = 3.0e17   # nm/sec
    input_photons = input_spec / (h * c / input_wv)
    objspec = np.interp(instrument.waves, input_wv, input_photons)
    src = objspec * magellan.Atel * texp * (instrument.waves) / instrument.R * \
        instrument.throughput * tel_throughput * atm      # point: e- from FULL source; extended: e-/arcsec2

    # Detector noise variance in ONE fibre (read + dark over the fibre's pixel profile).
    readnoise = instrument.sensor.rn
    dark = instrument.sensor.dark * texp
    pix_resel = (instrument.fiber.dFib) * (instrument.f_cam / instrument.f_col) / instrument.sensor.pixelsize
    det_var_perfib = (dark + readnoise ** 2) * np.ceil(pix_resel)

    # ---- extended source: surface brightness fills the fibre(s); no aperture loss ----
    if source == 'extended':
        n = max(int(nbin), 1)
        counts = n * src * instrument.fiber.Afib
        noise = np.sqrt(n * (sky_perfib + det_var_perfib))
        print("   extended source (surface brightness); co-adding {:d} fibre(s)".format(n))
        return counts, noise

    # ---- point source: seeing + fibre aperture (enclosed energy) ----
    if seeing is None:
        seeing = 0.8
        print("WARNING: no seeing given -- assuming {:.1f}\" FWHM. Pass seeing=<arcsec> for your "
              "conditions; the fibre aperture captures only part of a point source.".format(seeing))

    def _eval(N):
        r = aperture_radius(N, instrument.fiber.Afib)
        fenc = enclosed_energy(seeing, r, psf, beta)
        sig = src * fenc
        noise = np.sqrt(N * (sky_perfib + det_var_perfib))
        with np.errstate(invalid='ignore', divide='ignore'):
            snr = sig / noise
        return fenc, sig, noise, np.nanmedian(snr[np.isfinite(snr)])

    if aperture == 'single':
        N = 1
    elif isinstance(aperture, (int, np.integer)) and not isinstance(aperture, bool):
        N = max(int(aperture), 1)
    elif aperture == 'optimal':
        N = max(range(1, 20), key=lambda k: _eval(k)[3])   # aperture that maximises median band SNR
    else:
        raise ValueError("aperture must be 'single', 'optimal', or an integer number of fibres")

    fenc, counts, noise, _ = _eval(N)
    print("   point source: seeing {:.2f}\" ({}), aperture = {:d} fibre(s), enclosed energy = {:.2f}"
          .format(seeing, psf, N, fenc))
    return counts, noise
