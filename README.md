# llamas-etc
Public exposure time calculator for the Magellan LLAMAS Integral Field Spectrograph

Please review the supplied jupyter notebook LLAMAS_ETC_demo.ipynb for instructions on how to run the code.

## Graphical interface (optional)

If you prefer a GUI to the notebook, run the PyQt6 application:

```bash
pip install PyQt6        # one-time, if not already installed
python etc_gui.py
```

It's a single window wrapping the same `observe_spectrum` engine: set the input spectrum (defaults to the
bundled `SN1a_R20mag.fits`), exposure time, airmass, seeing, source type (point / extended surface
brightness), aperture, PSF, throughput model, and channels; the embedded plot shows counts and SNR vs
wavelength, and "Save results" writes a PNG + CSV.

## Throughput: as-measured (default) vs. theoretical

As of the first on-sky calibration, the ETC **defaults to the as-measured instrument throughput** derived
from spectrophotometric standard-star observations (the telescope+instrument response, with atmospheric
extinction removed). The pre-ship theoretical model (product of individual optical-element throughputs) is
still fully available:

```python
llamas_green = spec.Spectrograph('LLAMAS_GREEN')
llamas_green.build_model('llamas_green.def')                          # measured (default)
llamas_green.build_model('llamas_green.def', throughput_mode='theoretical')   # pre-ship model
```

The measured curves live in `COATINGS/measured_throughput_{blue,green,red}.txt` (2 columns: wavelength [nm],
throughput fraction). Where the measurement has no clean coverage — the channel edges and the broad red
telluric H2O complex (~890–990 nm) — the ETC automatically falls back to the theoretical curve. In measured
mode the telescope mirror reflectivity is already contained in the curve, so `observe.py` does not re-apply
it (this is handled automatically via `instrument.throughput_mode`). The as-measured throughput is currently
based on the may 2026 standards (GD108, Feige110); the absolute scale is provisional pending a review of
night-to-night transparency, but the wavelength shape is robust.

### Airmass (atmospheric extinction)

Both the measured and theoretical throughputs are referenced to **above the atmosphere**, so `observe_spectrum`
applies atmospheric extinction for your observation's airmass:

```python
counts, noise = observe.observe_spectrum(llamas_green, texp, wave_nm, flux)               # airmass 1.0 (warns)
counts, noise = observe.observe_spectrum(llamas_green, texp, wave_nm, flux, airmass=1.4)   # sec(z) = 1.4
```

If you don't pass `airmass` it **defaults to 1.0 (zenith) and prints a warning** — always pass your actual
sec(z). Extinction uses the Las Campanas curve in `COATINGS/lco_extinction.txt` and is applied to the source
as `10^(-0.4·k(λ)·airmass)` (strongest in the blue). The sky-background airglow is emitted high in the
atmosphere and is not extincted like a source, so it is left unmodified (its mild airmass dependence is a
future refinement).

### Point sources (seeing) and extended sources

For an **unresolved (point) source** (`source='point'`, the default; `input_spec` = flux density
erg/cm2/s/A), the ETC computes the fraction of light captured by the fibre aperture from the seeing PSF:

```python
counts, noise = observe.observe_spectrum(llamas_green, texp, wave_nm, flux,
                                         airmass=1.4, seeing=0.8, aperture='optimal')
```

`seeing` is the FWHM in arcsec (defaults to 0.8" with a warning); `aperture` is `'optimal'` (the fibre count
that maximises SNR — default), `'single'` (one fibre), or an integer number of fibres; `psf` is `'moffat'`
(default, index `beta`) or `'gaussian'`. Pass `seeing=0` to recover the pre-v1.1 "all light in one fibre"
behaviour.

For a **resolved (extended) source**, pass surface brightness (erg/cm2/s/A/arcsec²) with `source='extended'`;
each fibre samples SB × (fibre solid angle) with no aperture loss, and `nbin` co-adds fibres (SNR ∝ √nbin):

```python
counts, noise = observe.observe_spectrum(llamas_green, texp, wave_nm, sb, source='extended', nbin=1)
```

# Caveats

There are a few items to remember when interpreting results from the exposure time calculator:

1) **Point sources / seeing (now handled):** for an unresolved source the ETC computes the fraction of light captured by the fibre aperture from the seeing PSF, rather than assuming all light lands in one fibre. Set `seeing` (FWHM arcsec) and `aperture` (`'optimal'`/`'single'`/integer); pass `seeing=0` for the old full-capture behaviour. (See "Point sources (seeing) and extended sources" above.)

2) **Extended sources / surface brightness (now handled):** pass surface brightness (erg/cm2/s/A/arcsec²) with `source='extended'` and the ETC integrates over the fibre solid angle for you; use `nbin` to co-add fibres. Earlier versions required manually multiplying your surface-brightness spectrum by `fiber.Afib` and passing it as flux — that workaround is no longer necessary.

3) The ETC uses an average value for spectral R over the full instrument range and does not account for variations across the bandpass. This will be added in future releases.

4) The default throughput is now the **as-measured** on-sky telescope+instrument response from
spectro-photometric standard stars (may 2026: GD108, Feige110), replacing the pre-ship best-effort model as
the default (the theoretical model remains available via `throughput_mode='theoretical'`, and is used as an
automatic fallback where the standards give no clean measurement). The measured throughput runs ~0.6–0.75x
the pre-ship theoretical prediction, so exposure-time forecasts are now more realistic (and more
conservative) than earlier versions. Note the absolute scale is provisional pending a night-to-night
transparency review; the wavelength shape is robust.

If you develop new functionality for the ETC and would like to share this with future users, please contact the site administrator to request developer access to submit a pull request. We encourage regular users to share these improvements with the rest of the community!

Good luck with your observations and thanks for your interest in LLAMAS.
