#!/usr/bin/env python
"""Regenerate COATINGS/eso_moon_color.npz -- the ESO SkyCalc scattered-moonlight COLOUR template used by
observe.moon_sky_radiance() to colour the Krisciunas & Schaefer (1991) V-band moon level.

Requires network access to the ESO SkyCalc service and the skycalc_cli package (`pip install skycalc_cli`).
It queries the scattered-moonlight component (`flux_sml`, already photon radiance) at a set of moon-target
separations -- with hand-picked geometries that satisfy SkyCalc's |z - zmoon| <= rho <= z + zmoon
constraint -- normalises each spectrum to 1.0 at 550 nm, and stores color[nsep, nwave] plus the wavelength
and separation axes. Observatory '2400' (La Silla, ~2400 m) is the closest SkyCalc site to Las Campanas
(~2380 m). Separation dominates the moonlight colour (phase and airmass are second-order), so a
separation grid at full moon is used.

Run:  python tools/generate_moon_color.py
"""
import os
import tempfile

import numpy as np
from astropy.io import fits
from skycalc_cli.skycalc import SkyModel

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'COATINGS', 'eso_moon_color.npz')

# (moon-target separation [deg], moon altitude [deg], target airmass) -- each satisfies SkyCalc's geometry
GEOM = [(10, 80, 1.05), (30, 60, 1.15), (60, 60, 1.20), (90, 40, 1.40), (120, 20, 1.70), (150, 10, 3.00)]


def query_flux_sml(sep, alt, airmass):
    skm = SkyModel()
    skm.callwith({'airmass': airmass, 'moon_sun_sep': 180.0, 'moon_target_sep': float(sep),
                  'moon_alt': float(alt), 'observatory': '2400',
                  'wmin': 350.0, 'wmax': 1000.0, 'wdelta': 1.0, 'wgrid_mode': 'fixed_wavelength_step'})
    fn = tempfile.mktemp(suffix='.fits')
    skm.write(fn)
    t = fits.open(fn)[1].data
    os.remove(fn)
    return np.asarray(t['lam'], float), np.asarray(t['flux_sml'], float)


def main():
    wave, cols, seps = None, [], []
    for sep, alt, airmass in GEOM:
        lam, sml = query_flux_sml(sep, alt, airmass)
        if wave is None:
            wave = lam
        sml = np.interp(wave, lam, sml)
        cols.append(sml / sml[np.argmin(np.abs(wave - 550.0))])       # normalise colour at 550 nm
        seps.append(sep)
        print('  sep=%3d deg done' % sep)
    np.savez(OUT, wave_nm=wave, sep_deg=np.array(seps, float), color=np.array(cols))
    print('wrote %s  (color shape %s)' % (OUT, np.array(cols).shape))


if __name__ == '__main__':
    main()
