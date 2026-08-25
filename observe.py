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


def moon_sky_radiance(waves_nm, illum, sep_deg, alt_deg, airmass, k_V, T_sun=5777.0):
    """Scattered-moonlight sky background as an ADDED photon radiance [photons/m2/s/micron/arcsec2] on
    ``waves_nm`` -- the Krisciunas & Schaefer (1991) optical model. Zero for a new moon or a moon below
    the horizon.

        illum   : lunar illuminated fraction, 0 (new) .. 1 (full).
        sep_deg : moon-target angular separation [deg].
        alt_deg : moon altitude [deg].
        airmass : target airmass (sec z).
        k_V     : V-band atmospheric extinction [mag/airmass].

    KS91 gives the V-band sky-brightness increase; it is coloured here by a solar (``T_sun``) Planck
    spectrum with a Rayleigh (550nm/lambda)^4 blue boost, normalised at 550 nm -- so the moon background
    is bluer than the dark sky. This is the standard ETC-level approximation (the ESO Paranal sky model
    is more detailed: multiple scattering, aerosols)."""
    waves_nm = np.asarray(waves_nm, float)
    if illum is None or alt_deg is None or alt_deg <= 0.0:
        return np.zeros_like(waves_nm)
    illum = min(max(float(illum), 0.0), 1.0)
    alpha = np.degrees(np.arccos(2.0 * illum - 1.0))              # phase angle: 0=full, 180=new [deg]
    m_moon = -12.73 + 0.026 * alpha + 4.0e-9 * alpha ** 4         # KS91 moon V magnitude
    I_star = 10.0 ** (-0.4 * (m_moon + 16.57))                   # lunar illuminance
    f_rho = 10.0 ** 5.36 * (1.06 + np.cos(np.radians(sep_deg)) ** 2) + 10.0 ** (6.15 - sep_deg / 40.0)
    s = np.sin(np.radians(90.0 - alt_deg))
    X_moon = (1.0 - 0.96 * s * s) ** (-0.5)                      # KS91 airmass at the moon
    B_moon = f_rho * I_star * 10.0 ** (-0.4 * k_V * X_moon) * (1.0 - 10.0 ** (-0.4 * k_V * airmass))  # nanoLamberts
    if not np.isfinite(B_moon) or B_moon <= 0.0:
        return np.zeros_like(waves_nm)
    mu_V = (20.7233 - np.log(B_moon / 34.08)) / 0.92104          # V surface brightness [mag/arcsec2]
    # spectral shape: solar Planck (energy) x Rayleigh (550/lambda)^4, normalised at 550 nm
    h = 6.62607e-27; c = 2.99792e10; kB = 1.38065e-16
    def planck_energy(lam_nm):
        lam_cm = np.asarray(lam_nm, float) * 1.0e-7
        return lam_cm ** -5 / (np.exp(h * c / (lam_cm * kB * T_sun)) - 1.0)
    shape = (planck_energy(waves_nm) / planck_energy(550.0)) * (550.0 / waves_nm) ** 4
    f_lambda = 3.631e-9 * 10.0 ** (-0.4 * mu_V) * shape          # erg/s/cm2/A/arcsec2 (V zeropoint x shape)
    photons_cm2_A = f_lambda * (waves_nm * 1.0e-7) / (h * c)     # photons/s/cm2/A/arcsec2
    return photons_cm2_A * 1.0e8                                 # -> photons/m2/s/micron/arcsec2


def observe_spectrum(instrument, texp, input_wv, input_spec, airmass=None,
                     source='point', seeing=None, aperture='optimal', psf='moffat', beta=2.5,
                     nbin=1, moon_illum=None, moon_sep=90.0, moon_alt=45.0,
                     skyfile="eso_newmoon_radiance.txt", extfile="lco_extinction.txt"):
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

    Moonlight (optional): pass ``moon_illum`` (illuminated fraction 0..1) to add a scattered-moonlight sky
        background (Krisciunas & Schaefer 1991) on top of the dark sky. ``moon_sep`` = moon-target
        separation [deg], ``moon_alt`` = moon altitude [deg]. Default (moon_illum=None) is the dark sky.
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

    # Optional scattered-moonlight background (Krisciunas & Schaefer 1991), added to the dark sky.
    if moon_illum is not None:
        try:
            k_V = float(np.interp(550.0, ke['wave_nm'], ke['k']))
        except Exception:
            k_V = 0.12
        moon = moon_sky_radiance(instrument.waves, moon_illum, moon_sep, moon_alt, airmass, k_V)
        sky = sky + moon
        tag = ('below horizon/new -> no added background' if not np.any(moon > 0)
               else 'added (bluer than dark sky)')
        print("   moon: illum={:.2f}, sep={:.0f} deg, alt={:.0f} deg -> {}"
              .format(min(max(moon_illum, 0.0), 1.0), moon_sep, moon_alt, tag))

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
