%
% Contents file for package: AstroSpec
% Created: 29-Dec-2015
%---------
% accretion_disk.m :  Calculate the theoretical spectrum of a optically thick, thin accretion disk given the mass of the central object, disk radius, and the accretion rate.
% accretion_disk_mag_c.m :  Calculate the magnitude, in a given filter, of a optically-thick thin accretion disk model.
% add_specphot_stand.m :  Add a spectrophotometric standard to list of standard stars.
% band_spectrum.m :  Calculate an un-normalized Band-spectrum (Band et al. 1993).
% black_body.m :  Calculate a black-body spectrum given its temperature.
% blackbody_bolmag.m :  Calculate the bolometric magnitude of a black body spectrum, given its temperature, radius and distance.
% blackbody_flux.m :  Calculate the flux, in a given wavelength range, of a black-body, given its temperature, radius and distance.
% blackbody_mag_c.m :  Calculate the magnitude, in a given filter, of a black-body, given its temperature, radius and distance.
% brightness_temp.m :  Calculate the brightness temperature.
% calibrate_spec_using_phot.m :  Given a spectrum and a light curve in one or more bands, applay synthetic photometry to the spectrum to calibrate it against the light curves. Optionally, if RA, Dec are provided, dereddned the spectrum for Galactic extinction. Furthermore, if redshift is provided then the spectrum will be converted to the rest frame, and if additional host extinction is given then another deredenning at the rest frame wavelength will be applied. If more than one filter is provided than the spectrum calibration will be performed using a polynomial fitting when the deg of the polynom is N-1, where N is the good filters.
% ccf_spectra.m :  Cross-Correlate two spectra. For each redshift, the correlation is calculated by applying a redshift to the second spectrum (Template), and correlate it with the first spectrum.
% chi2_spectra.m :  Perform a chi^2 fitting between a spectrum and a template spectrum as a function of redshift.
% conv1_gaussian.m :  Convolve a 1-D signal with a Gaussian kernel.
% convert_flux.m :  Convert between different flux units
% eq_temp.m :  Calculate the eqilibrium temperature of a body illuminated by a black-body radiation.
% extinction.m :  Given the E_{B-V} and the wavelength or filter name calculate the extinction in magnitude. The program works in the 0.1-2 micron range. The program is using the Cardelli, Clayton, Mathis (1989) model. See also old version: optical_extinction.m.
% find_shift_scale_spec.m :  Given a spectrum and a reference spectrum [Wavelength, Flux] find the shift (a constant added to the data) and/or a scaling (multiply the data by a constant) that relats the data and the reference by minimizing the rms. Note that the shift is applaied first. Alternatively, the program can find a factor as function of wavelength that such: Ref*Factor = Spec. This is done by fitting a polynomial to the log of the ratio between the spectra and reference.
% fit_blackbody.m :  Fit set of magnitudes to a black-body spectrum, and derive its best fit parameters.
% fit_blackbody_spec.m :  Fit a spectrum to a black-body spectra.
% fit_bolometric_flux.m :  Estimate the Total bolometric flux according to given magnitudes of different filters.
% fit_extinction.m :  Given a spectrum and a template, find the best extinction curve that should be applied to the spectrum in order to best fit the template.
% fit_photoz.m :  Photometric redshift \chi^2 estimator. Given a photometric measurment of a source and a spectrum, calculate the \chi^2 fit between the spectrum as a function of redshift and the photometric measurements.
% fit_specline.m :  
% fun_gauss.m :  Calculating a Gaussian function of the form: Y = Amplitude*exp( (X-W0)^2/(2*Sigma^2) ) Optionaly, convolve the the result with a Gaussian. This function can be used by fitting functions like nlinfit_my.m and fit_specline.m
% fun_lorentzian.m :  Calculating a Lorentzian function of the form: Y = D.*Gamma./(pi.*( (X-X0).^2 + Gamma.^2 )) Optionaly, convolve the the result with a Gaussian. This function can be used by fitting functions like nlinfit_my.m and fit_specline.m
% get_filter.m :  Search and get astronomical Filter information and transmission curve.
% get_gaia_synspec.m :  get a synthetic stellar spectrum from the local GAIA spectral library. Spectra are in the range 2500-10500A and 1A resolution. Assuming alpha enhanement 0, and micro-turbulence 2km/s.
% get_lines.m :  Search spectral lines by wavelength. (for search by name see: get_lines1.m) The list of lines contains 46663 spectral lines (Reader et al. 1980; 1981) for 99 atomic species. Neutral through quadruply ionized atoms are tabulated. 
% get_lines1.m :  Search for a spectral line, by name or by eavelength, among a list of selected spectral lines.
% get_spectra.m :  Get a template spectrum from spectra library, or alternatively a use supplied spectrum. Optionally, apply redshift, extinction, atmospheric extinction, and Telluric absorptions.
% hydrogen_lines.m :  Calculate the vacum wavelength of Hydrogen lines, given their shell numbers.
% interp_mag.m :  Given magnitude of an object find the best fit spectra (from a library of templates), and calculate the magnitude of the spectra in additional bands.
% ionization_potential.m :  Returm the ionization potential for a given element and ionization level.
% kcorr.m :  Calculate k-correction. Given a spectrum two filters, and their redshifts calculate the k-correction of the first filter minus the second filter. This is calculated by shifting the spectrum to z1, measuring the synthetic magnitude in the first filter, then shiting to z2 and measuring the synthetic magnitude in the second filter.
% luptitude.m :  Convert flux to luptitudes (asinh magnitudes).
% matchspec.m :  A GUI utility to inspect and match spectrum with templates.
% optical_extinction.m :  Given the E_{\lambda_1-\lambda_2} (e.g., E_{B-V}) calculate the extinction in magnitude A_{\lambda_3}. The program works in the 0.1-2 micron range. The program is using the Cardelli, Clayton, Mathis (1989) or Allen models. \lambad can be specified using filters names: Johnson: 'U','B','V','R','I','J','H','K' POSS : 'O','E' SDSS : 'u','g','r','i','z'
% scale_spectrum.m :  Scale spectrum by shift and stretch or a wavelength dependent factor (See also: find_shift_scale_spec.m).
% search_specphot_stand.m :  Search a spectroscopic standard star, by coordinates or by name in the SpecPhot_Stand.mat structure. The search by name option is case insensitive and blanks insensitive. Moreover, the name search can be either exact or by substring.
% shift2vel.m :  Calculate the velocity from the red/blue shift (z).
% shift_spec.m :  Transform a spectrum from the observed frame to the rest frame. If the redshift is negative then transform from the rest frame to the observed frame.
% sky_ebv.m :  Calculate the extinction from a local copy of the Schlegel, Finkbeiner & Davis (1998) extinction maps.
% spec_photon_counts.m :  Given a spectrum and the effective area of an instrument as a function of wavelength, calculate the the total recieved flux and the photons count rate in the instrument.
% star_sptype_color.m :  Given a star spectral type and luminosity class, get the star color between any two filters.
% synphot.m :  Calculate synthetic photometry of a spectrum.
% vel2shift.m :  Calculate the red/blue shift (z) from velocity.
% wein.m :  Apply Wein law - return the peak wavelength of a black body at a given temperature.
% wget_gaia_synspec.m :  wget a synthetic stellar spectrum from the GAIA spectral library, in the range 2500-10500A and 1A resolution. Assuming alpha enhanement 0, and micro-turbulence 2km/s.
% xray_abs.m :  Given the neutral Hydrogen column density, calculate the bound-free attenuation of X-rays as a function of wavelength in the ISM. The program assumes abundences from Ebihara (1982). Absorption is due to neutral species only. Adopted from Zombeck (1990).
% zodiac_spectrum.m :  Return the zodiac spectrum as adopted from the HST STIS handbook. The high zodiacal ligh is defined where V=22.1 mag/arcsec^-2.
