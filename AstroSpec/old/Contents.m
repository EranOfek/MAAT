%
% Contents file for package:  AstroSpec 
% Created on 09-Apr-2007 
%-----------------------------------------------------
%
% accretion_disk - Calculate an accretion disk spectrum given the mass of
%                 the central object, disk radii, and the accretion rate.
% band_spectrum - Calculate an un-normalized Band-spectrum (Band et al. 1993).
% black_body    - Calculate a black-body spectrum given its temperature.
% blackbody_mag - Calculate the magnitude, in a given bad, of a black body
%                 given its temperature, radius and distance.
%                 Work only in the range 100-10^6 K, for other temperatures
%                 and additional filters use: blackbody_mag_c.m
% blackbody_mag_c - Calculate the magnitude, in a given filter, of a black-body,
%                 given its radius and distance.
% calc_synflux  - Calculate synthetic flux in a given filter from
%                 a spectrum. The results are non zero-point calibrated.
% ccf_spectra   - Cross-Correlate two spectra. For each redshift, the
%                 correlation is calculated by applying a redshift to the
%                 second spectrum (Template), and correlate it with the
%                 first spectrum.
% chi2_blackbody_mag - Given a list of observed magnitudes in several bands,
%                 perform a \chi^2 fit with a black-body spectrum with a
%                 given temperatures and extinctions.
% fit_photoz    - Photometric redshift \chi^2 estimator. Given a photometric
%                 measurment of a source and a spectrum, calculate the \chi^2
%                 fit between the spectrum as a function of redshift and the
%                 photometric measurements.
% get_filter    - Search and get astronomical Filter information
%                 and transmission curve.
% get_lines     - Search spectral lines by wavelength.
%                 (for search by name see: get_lines1.m)
%                 The list of lines contains 46663 spectral lines
%                 (Reader et al. 1980; 1981) for 99 atomic species.
%                 Neutral through quadruply ionized atoms are tabulated. 
% get_lines1    - Search for a spectral line, by name or by eavelength,
%                 among a list of selected spectral lines.
% get_spectra   - Get a template spectrum from spectra library,
%                 or alternatively a use supplied spectrum.
%                 Optionally, apply redshift, extinction, atmospheric extinction,
%                 and Telluric absorptions.
% interp_mag    - Given magnitude of an object find the best fit spectra
%                 (from a library of templates), and calculate the
%                 magnitude of the spectra in additional bands.
% mag_jy_conv   - Calculate conversion factor between magnitude, flux and
%                 specific flux units (e.g., mJy).
%                 Allows for user supplied filter and spectral slope.
% matchspec     - A GUI utility to inspect and match spectrum with templates.
% optical_extinction - Given the E_{\lambda_1-\lambda_2} (e.g., E_{B-V})
%                 calculate the extinction in magnitude A_{\lambda_3}.
%                 The program works in the 0.1-2 micron range.
%                 The program is using the Cardelli, Clayton, Mathis (1989)
%                 or Allen models.
%                 \lambad can be specified using filters names:
%                 Johnson: 'U','B','V','R','I','J','H','K'
%                 POSS   : 'O','E'
%                 SDSS   : 'u','g','r','i','z'
% plot_sdss_tc  - Plot SDSS filters transmission curves
% scale_spectrum - Scale spectrum by shift and stretch or wavelength
%                 dependent factor (See also: find_shift_scale_spec.m).
% shift2vel     - Calculate the velocity from the red/blue shift (z).
% shift_spec    - Transform a spectrum from the observed frame to
%                 the rest frame.
% sky_ebv       - Calculate the extinction from a local copy of the 
%                 Schlegel, Finkbeiner & Davis (1998) maps.
% vel2shift     - Calculate the red/blue shift (z) from velocity.
% wein          - Apply Wein law - return the peak wavelength of a
%                 black body at a given temperature.
