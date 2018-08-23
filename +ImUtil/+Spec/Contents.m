%
% Contents file for package: ImSpec
% Created: 29-Dec-2015
%---------
% atmospheric_ext.m :  Correct a spectrum for atmospheric extinction.
% def_bitmask_specpipeline.m :  The spectroscopic pipeline bit mask definition. Given the Bit mask name return the bit mask index.
% def_instrument_arc.m :  Definitions of how to identify spectroscopic arc accoring to the image header information. Given an instrument, this function return a header keyword name that contains the arc information, as well as the values of this keyword if the image is an arc image, and the translation of these various values to standard arc name.
% eq_width.m :  Given a spectrum calculate the equivalent width, flux, and center of selected spectral lines. The script has interactive and noninteractive modes. * Interactive mode: The user select (left click) pairs of wavelength positions marking the continuum around each line (right click to abort). * Noninteractive mode: If a second input argument is given then the lines properties are calculated noninteractively. This function is obsolete, use: fit_specline.m instead.
% find_contrast_peaks.m :  Given a list of [X, Y], find the position of local maxima in the list which are above a certain contrast above the local rms. The rms is defined as the 68 percentile of the extramum Y-value distribution, while the contrast is defined as the Y-value of the local maxima subtract from the preceeding local minima.
% find_src1d.m :  Find sources (e.g., local maxima above the noise) in 1-D data. Find their location, S/N, width, and local background.
% is_arc_image.m :  Given a list of FITS images or SIM, look for arc (wavelength calibration) images. The search is done by looking for specific keyword values in the image headers.
% is_stdstar_image.m :  Given a list of FITS images or SIM, look for standard stars spectra. The search is done by comparing the RA/Dec or/and object name keywords with the database of spectroscopic standard stars.
% read_fits_spec.m :  Get a spectrum from a FITS file. Currently support on CTYPE1='linear'.
% spec_classify_images.m :  Classify a list of images given their header keywords. For each imagee the program retrieve or calculate the value of some header keywords or some constants. Moreover, for another set of keywords the retrival mode may depends on the value of a specific keyword (e.g., specifying if a blue or red arm image).
% spec_coadd1d.m :  
% spec_collapse_dispaxis.m :  Given an image/s collapse each image to 1 D vector along one of the dimensions (e.g., "dispersion direction").
% spec_convolve_resolution.m :  Given a spectrum and a constant resolution term (dLambda/Lambda), convolve a filter which width equal to the resolution, with the spectrum. This function can be used to degrade a spectrum to a specific resolution.
% spec_dbsp_bin.m :  
% spec_extract_2d.m :  Given an image and a trace, cut a sub image along the trace position in which the trace is aligned along the X-axis.
% spec_fit1d_psf.m :  Given a traced and sky subtracted spectrum, extract the intensity of a single spectrum, along the dispersion axis, by fitting a 1-dimensional numerical PSF. The PSF is constructed by binning the spectrum.
% spec_fluxcalib.m :  Given the observed spectrum of a flux standard, and the flux standard calibrated spectrum outside the atmosphere, attempt to find the transmission function of the spectrograph and telescope.
% spec_get_arc.m :  Get a spectroscopic arc (template and lines list) from the SpecArcs.mat database, given the arc name.
% spec_response.m :  Given the observed spectrum of a flux standard, and the flux standard calibrated spectrum outside the atmosphere, attempt to find the transmission function of the spectrograph and telescope.
% spec_skysub.m :  Subtract the sky from a 2-dimensional spectrum.
% spec_skysub1.m :  Subtract the sky from a 2-dimensional spectrum.
% spec_stitch.m :  Given two spectra of a single object taken at different wavelength (e.g., blue side and red side), with or without an overlap between the spectra, this program find an optimal flux ratio shift between the two spectra, and return the stitched spectrum.
% spec_telluric_template.m :  Generate a Telluric template. Give a spectrum of the ratio between the calibrated observed spectrum and the standard star theoretical spectrum, and a vector of flags indicating if a wavelength is a Telluric line, return a template of the Telluric lines.
% spec_trace.m :  Trace and extract spectrum from a 2-dimensional image.
% spec_trace_extract_wave.m :  A wrapper around spec_trace.m, spec_extract 2/1-d and automatic wavelength calibration.
% spec_trace_int.m :  Trace a spectrum in a single image by tracing the peak of the intensity pixel by pixel in the dispersion direction.
% spec_wavecalib_lines.m :  Given a spectrum with an approximate wavelength solution, identify lines in the spectrum and attempt to match them with list of spectral lines in template spectra (e.g., Arcs or sky lines). Then, refine the wavelength calibration by fitting a function (e.g., polynomial) to the lines in the spectrum to match the lines in the template.
% spec_wavecalib_xcorr.m :  Perform automatic wavelength calibration for 1-D spectrum, by cross corelating it with a template spectrum. The cross-correlation is done using various trial scales, so the program optimize for both the scale and shift. After the cross-correlation, refine the wavelength calibration by line matching.
