%
% Contents file for package: SEDM
% Created: 29-Dec-2015
%---------
% sedm_abswavecalib.m :  Find an absolute wavelength calibration solution, relative to a calibrated template, for selected segments. The solutions are found by searching for a distortion transformation that match between the template and the observed spectra. This is later used by sedm_wavecalib.m to solve for the wavelength calibration for all the segments.
% sedm_abswavecalib1.m :  Find an absolute wavelength calibration solution, relative to a calibrated template, for selected segments. The solutions are found by searching for a distortion transformation that match between the template and the observed spectra. This is later used by sedm_wavecalib.m to solve for the wavelength calibration for all the segments. This is an older version of sedm_abswavecalib.m that may be usefull sometimes.
% sedm_abswavecalib_indiv.m :  Find an absolute wavelength calibration solution, relative to a calibrated template, for selected segments. The solutions are found by searching for a distortion transformation that match between the template and the observed spectra. This is later used by sedm_wavecalib.m to solve for the wavelength calibration for all the segments.
% sedm_bias_subtraction.m :  Subtract bias from SEDM image/s.
% sedm_copy_wavecalib.m :  Copy wavelength calibration information from a SegmentsInfo structure array of an Arc (i.e., after sedm_wavecalib.m) to that of a science SegmentsInfo structure array.
% sedm_extract_spec_cube.m :  Extract a source and background spectrum from a SEDM cube.
% sedm_extract_spexcell.m :  Given a SEDM SegmentsInfo structure returned by sedm_generate_spexcell_segmentation.m and a science or calibration image, search for the exact position and slope of the trace in each spexcell using the SegmentsInfo as a first guess.
% sedm_generate_spexcell_segmentation.m :  Use SEDM dome flat images to generatre segmentation map of the individual spexcells (i.e., each segment is marked by a serial number).
% sedm_plot_seg_int.m :  Calculate for each spexcell its median value within a givenb wavelength range and store it in SegmentsInfo, and plot the median intensity of each spexcell as a function of its mean X/Y position.
% sedm_plottrace.m :  plot a spexcell trace/s on image.
% sedm_prep_nightlog.m :  Given a list of SEDM images, prepare a log file of images with basic information derived from the image headers.
% sedm_refine_slopepos.m :  Given a SEDM SegmentsInfo structure returned by sedm_generate_spexcell_segmentation.m and a science or calibration image, search for the exact position and slope of the trace in eac spexcell using the SegmentsInfo as a first guess.
% sedm_segments2cube.m :  Given a SegmentsInfo structure that contains in each element the wavelength calibrated spectrum of each spexcell, generate a cube of flux values as a function of X position, Y position and wavelength.
% sedm_wavecalib.m :  Find and apply the wavelength calibration to all the spexcells. The function ...
% sedm_wavelength_xcorr.m :  Given two spectra and a distortion transformation describing the relation between the wavelength of the two spectra, return the negative of the cross-correlation of the two spectra after applying the transformation.
