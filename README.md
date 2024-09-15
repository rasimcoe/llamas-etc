# llamas-etc
Public exposure time calculator for the Magellan LLAMAS Integral Field Spectrograph

Please review the supplied jupyter notebook LLAMAS_ETC_demo.ipynb for instructions on how to run the code.

# Caveats

There are a few items to remember when interpreting results from the exposure time calculator:

1) The ETC assumes that all of the light goes down a single fiber, i.e. that the source is unresolved in the 0.75" spaxel. If your source is resolved, or the seeing is lousy, you will need to adjust these outputs accordingly, and we have not yet implemented that functionality (though it would be straightforward to split the light across N fibers/spaxels).  

2) The ETC uses units of flux, and not surface brightness, but surface brightness is the more appropriate unit for resolved sources.  Future ETC versions may include surface brightness calculations, but intrepid users can implement this on their own using the following hints. Surface brightness can often be expressed in units of erg/cm2/s/A/square arcsec, or can be converted into these units from magnitudes per square arcsec or your unit of choice. It is your responsibility to convert into the first set of units (erg/cm2/s/A/sq.arcsec). Then, the subtended area of a fiber is provided as an attribute in the spectrograph object: for example, llamas_blue.fiber.Afib (in square arcseconds). Multiply your surface brightness spectrum by this value, and then input that to the observe.observe_spectrum subroutine to output a correct SNR calculation!

3) The ETC uses an average value for spectral R over the full instrument range and does not account for variations across the bandpass. This will be added in future releases.

4) Once again, we remind users that this is a best-effort estimate based on multiplying the measured and model curves for optical glasses, coatings, mirrors, the fibers, and geometric shadows.  Users are urged to plan conservatively for early observing runs, until the team has an opportunity to observe spectro-photometric standard stars on the sky. When those observations are taken we will update this repository with as-measured throughputs so that future observers can forecast with confidence.

If you develop new functionality for the ETC and would like to share this with future users, please contact the site administrator to request developer access to submit a pull request. We encourage regular users to share these improvements with the rest of the community!

Good luck with your observations and thanks for your interest in LLAMAS.
