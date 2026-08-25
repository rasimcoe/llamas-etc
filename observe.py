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


_MOON_COLOR = None


def _moon_color_shape(waves_nm, sep_deg):
    """ESO SkyCalc scattered-moonlight colour (``flux_sml`` normalised to 1.0 at 550 nm) interpolated to
    the requested moon-target separation and wavelengths, from the bundled grid
    COATINGS/eso_moon_color.npz. Falls back to a flat (grey) spectrum if the grid is missing."""
    global _MOON_COLOR
    if _MOON_COLOR is None:
        p = 'COATINGS/eso_moon_color.npz'
        p = p if os.path.isfile(p) else os.environ.get('COATINGS_PATH', './COATINGS/') + 'eso_moon_color.npz'
        try:
            d = np.load(p)
            _MOON_COLOR = (np.asarray(d['wave_nm'], float), np.asarray(d['sep_deg'], float),
                           np.asarray(d['color'], float))
        except Exception:
            _MOON_COLOR = False
    if not _MOON_COLOR:
        return np.ones_like(np.asarray(waves_nm, float))
    gw, gs, gc = _MOON_COLOR
    s = float(np.clip(sep_deg, gs.min(), gs.max()))
    j = int(np.clip(np.searchsorted(gs, s) - 1, 0, len(gs) - 2))
    w = (s - gs[j]) / (gs[j + 1] - gs[j]) if gs[j + 1] != gs[j] else 0.0
    col = (1.0 - w) * gc[j] + w * gc[j + 1]                       # colour at this separation (grid wave grid)
    return np.interp(waves_nm, gw, col, left=col[0], right=col[-1])


def moon_sky_radiance(waves_nm, illum, sep_deg, alt_deg, airmass, k_V):
    """Scattered-moonlight sky background as an ADDED photon radiance [photons/m2/s/micron/arcsec2] on
    ``waves_nm``. Zero for a new moon or a moon below the horizon.

        illum   : lunar illuminated fraction, 0 (new) .. 1 (full).
        sep_deg : moon-target angular separation [deg].
        alt_deg : moon altitude [deg].
        airmass : target airmass (sec z).
        k_V     : V-band atmospheric extinction [mag/airmass].

    The V-band brightness LEVEL follows Krisciunas & Schaefer (1991); the spectral COLOUR is taken from
    the ESO SkyCalc sky model (Jones et al. 2013 scattered moonlight; bundled grid
    COATINGS/eso_moon_color.npz), interpolated by moon-target separation and normalised at 550 nm. The
    ESO colour is far less steeply blue than a single-scattering Rayleigh law -- important for the LLAMAS
    blue channel. (Separation dominates the colour; phase and airmass are second-order.)"""
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
    mu_V = (20.7233 - np.log(B_moon / 34.08)) / 0.92104          # V surface brightness [mag/arcsec2] (KS91)
    shape = _moon_color_shape(waves_nm, sep_deg)                 # ESO SkyCalc PHOTON colour, = 1.0 at 550 nm
    h = 6.62607e-27; c = 2.99792e10                              # erg s ; cm/s
    # KS91 level -> V PHOTON radiance in the ETC's sky units [photons/m2/s/micron/arcsec2] (f_lambda V
    # zeropoint / photon energy at 550 nm), then apply the ESO photon colour. flux_sml is already photon
    # radiance, so there is NO further f_lambda->photon (x lambda) conversion.
    R_V = 3.631e-9 * 10.0 ** (-0.4 * mu_V) * (550.0e-7) / (h * c) * 1.0e8
    return R_V * shape


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
