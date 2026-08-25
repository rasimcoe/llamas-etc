"""Basic sanity/regression tests for the LLAMAS ETC v1.1 features.

Run from the repository root (the modules import each other by bare name and read COATINGS/):
    python tests/test_observe.py        # or:  pytest -q
"""
import os
import sys

import numpy as np

# make the repo root importable + resolvable regardless of where pytest is invoked
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)
os.environ.setdefault('COATINGS_PATH', './COATINGS/')

import spectrograph as spec          # noqa: E402
import observe                       # noqa: E402
from astropy.table import Table      # noqa: E402


def _setup():
    g = spec.Spectrograph('LLAMAS_GREEN')
    g.build_model('llamas_green.def')
    sn = Table.read('SN1a_R20mag.fits')
    return g, np.array(sn[sn.colnames[0]]), np.array(sn[sn.colnames[1]])


def _med(a):
    return float(np.nanmedian(a[np.isfinite(a)]))


def test_enclosed_energy_limits():
    # perfect seeing -> everything enclosed; grows with radius; -> 1 at large radius
    assert observe.enclosed_energy(0.0, 0.4) == 1.0
    assert observe.enclosed_energy(0.8, 0.4) < observe.enclosed_energy(0.8, 1.2)
    assert observe.enclosed_energy(0.8, 20.0) > 0.99


def test_seeing_zero_recovers_full_capture():
    g, wv, fl = _setup()
    c_full, _ = observe.observe_spectrum(g, 1200, wv, fl, airmass=1.2, seeing=0.0)
    c_seeing, _ = observe.observe_spectrum(g, 1200, wv, fl, airmass=1.2, seeing=0.8, aperture='single')
    assert _med(c_full) > _med(c_seeing)                 # aperture loss reduces captured light


def test_extended_matches_point_full_capture():
    # surface brightness SB = flux/Afib must equal a fully-captured point source of that flux
    g, wv, fl = _setup()
    c_full, _ = observe.observe_spectrum(g, 1200, wv, fl, airmass=1.2, seeing=0.0)
    c_ext, _ = observe.observe_spectrum(g, 1200, wv, fl / g.fiber.Afib, airmass=1.2, source='extended')
    assert abs(_med(c_ext) / _med(c_full) - 1.0) < 1e-6


def test_optimal_beats_single():
    g, wv, fl = _setup()
    c_o, n_o = observe.observe_spectrum(g, 1200, wv, fl, airmass=1.2, seeing=0.9, aperture='optimal')
    c_s, n_s = observe.observe_spectrum(g, 1200, wv, fl, airmass=1.2, seeing=0.9, aperture='single')
    assert _med(c_o / n_o) >= _med(c_s / n_s)            # optimal aperture is at least as good


def test_higher_airmass_lower_counts():
    g, wv, fl = _setup()
    c1, _ = observe.observe_spectrum(g, 1200, wv, fl, airmass=1.0, seeing=0.8)
    c2, _ = observe.observe_spectrum(g, 1200, wv, fl, airmass=2.0, seeing=0.8)
    assert _med(c2) < _med(c1)


def test_full_moon_raises_sky_noise():
    g, wv, fl = _setup()
    _, n_dark = observe.observe_spectrum(g, 1200, wv, fl, airmass=1.2, seeing=0.8)
    _, n_moon = observe.observe_spectrum(g, 1200, wv, fl, airmass=1.2, seeing=0.8,
                                         moon_illum=1.0, moon_sep=60, moon_alt=60)
    assert _med(n_moon) > _med(n_dark)                   # bright moon brightens the sky


def test_moon_below_horizon_equals_dark():
    g, wv, fl = _setup()
    _, n_dark = observe.observe_spectrum(g, 1200, wv, fl, airmass=1.2, seeing=0.8)
    _, n_down = observe.observe_spectrum(g, 1200, wv, fl, airmass=1.2, seeing=0.8,
                                         moon_illum=1.0, moon_sep=60, moon_alt=-5)
    assert np.allclose(n_down, n_dark)                   # moon below the horizon adds nothing


def test_moon_background_is_bluer():
    m = observe.moon_sky_radiance(np.array([400.0, 900.0]), 1.0, 60, 60, 1.2, 0.12)
    assert m[0] > m[1]                                   # more moon background in the blue


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for fn in fns:
        fn()
        print("PASS", fn.__name__)
    print("all %d tests passed" % len(fns))
