# llamas-etc
Public exposure time calculator for the Magellan LLAMAS Integral Field Spectrograph

Please review the supplied jupyter notebook LLAMAS_ETC_demo.ipynb for instructions on how to run the code.

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

# Caveats

There are a few items to remember when interpreting results from the exposure time calculator:

1) The ETC assumes that all of the light goes down a single fiber, i.e. that the source is unresolved in the 0.75" spaxel. If your source is resolved, or the seeing is lousy, you will need to adjust these outputs accordingly, and we have not yet implemented that functionality (though it would be straightforward to split the light across N fibers/spaxels).  

2) The ETC uses units of flux, and not surface brightness, but surface brightness is the more appropriate unit for resolved sources.  Future ETC versions may include surface brightness calculations, but intrepid users can implement this on their own using the following hints. Surface brightness can often be expressed in units of erg/cm2/s/A/square arcsec, or can be converted into these units from magnitudes per square arcsec or your unit of choice. It is your responsibility to convert into the first set of units (erg/cm2/s/A/sq.arcsec). Then, the subtended area of a fiber is provided as an attribute in the spectrograph object: for example, llamas_blue.fiber.Afib (in square arcseconds). Multiply your surface brightness spectrum by this value, and then input that to the observe.observe_spectrum subroutine to output a correct SNR calculation!

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
