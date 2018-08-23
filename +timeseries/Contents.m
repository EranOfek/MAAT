%
% Contents file for package: timeseries
% Created: 29-Dec-2015
%---------
% arp.m :  Calculate the autoregessive process of order p. Moddeling an evenly spaced time series, z(t), with: z(t) = sum_i[a(i)*z(t-i)] + err(t) 
% bin_by_eye.m :  Bin light curve by eye. The function plots the light curve and let the user to define bins using the mouse. The user mark (with the mouse) the beginning and end points of each bin. The left and right limits of each user defined bin are marked on the plot using cyan and red dashed lines, respectively.
% binning.m :  Binning a timeseries using equally spaced bins. In each bin, calculate the mean, median, std, skewness of the data.
% ccf_diff.m :  Given two equally sapced time series of the same length, calculate the difference between points of distance N and calculate the correlation between this differences in the two serieses. Return the correlation as a function of distance N. This is usefull if trying to estimate if the correlation between two series is due to high or low frequency. Supports partial correlations.
% ccf_fft.m :  Calculate the normalized cross correlation function of two evenly spaced timeseries using fft.
% ccf_o.m :  Cross correlation function for two sets of equally spaced time series. 
% cusum.m :  cumulative sum (CUSUM) chart to detect non stationarity in the mean of an evenly spaced series.
% dcf.m :  Discrete cross-correlation function and structure function for two sets of unequaly spaced stationary time series.
% extract_phase.m :  Given a time series, extract observation made in a given phases range. The phase and epoch is defined by the user.
% folding.m :  Folding a time series into a period, and calculate the phase for each data point in the time series.
% folding_solarsys.m :  Folding a time series of a solar system object into a period, and calculate the phase for each data point in the time series. Taking into account ligh travel effect, phase correction due to phase angle, and phase angle changes in brightness.
% hjd.m :  Convert Julian Day (UTC) to Helicentric/Barycentric Julian Day for geocentric observer.
% matched_filter.m :  Run a matched filter (North filter) for a 1-D series and return the possible peaks (signal) above a given threshold.
% period.m :  Periodicy search in a time series. This function can be used to run different types of periodicity searches and window functions, and automatically select the range of frequencies to search.
% period_fft.m :  Calculate power spectrum for evenly spaced time series using fast Fourier transform. See also period.m
% period_norm.m :  Calculate the normalized normal power spectrum of a times series. See period.m for a more flexiable function.
% period_norm_order2.m :  Calculate the normalized power spectrum of the second order of a times series (i.e., assuming that the first and second harmonies have the same amplitude). See period.m for a more flexiable function.
% period_norm_solarsys.m :  Calculate the un-normalized normal power spectrum of a times series for a solar system object orbiting the Sun. Taking into account ligh travel effect, phase correction due to phase angle, and phase angle changes in brightness.
% period_normnl.m :  Calculate the classical (Lomb) power spectrum of a time series using no loops. See period.m for a more flexiable function.
% period_scargle.m :  Calculate the un-normalized Scargle power spectrum of a times series. See period.m for a more flexiable function.
% periodia.m :  Classical power-spectrum of a non-evenly spaced time series. The power spectrum is normalized by the variance of the data. THIS FUNCTION IS BEING REPLACED BY period.m
% periodis.m :  Lomb-Scargle periodigram. Calculate the Lomb-Scargle power spectrum for a time series. THIS FUNCTION IS BEING REPLACED BY period.m
% periodmulti_norm.m :  Calculate the normal power spectrum of a set of multiple times series which have common times. This program run the power spectrum simultaneously for all time serieses and its is faster than running the power spectrum on a single time series at a time.
% perioent.m :  Periodogram using information entropy. For each trail period, the phase-magnitude space is divided into m by m cells, and the information entropy is calculated. Then the probability to find observations in each square is MU(i) and the Entropy is (-Sigma{MU*ln(MU)}). The output Entropy is normalized by ln(m*m).
% plot_period_folder.m :  Given the power spectrum and light curve, let the user to interactively select peaks in the power spectrum, and present the folded (and binned) light curve for the selected periods.
% polysubs.m :  Subtract polynomial from a data set (no errors).
% runmean1.m :  running mean on a 1-D vector.
% sf_interp.m :  Interpolate a time series, and estimate the errors due to the interpolation using the structure function of the time series. The error in each interpolated point is calculated by adding in quadrature the error of the neastest point with the amplitude of the stracture function at the the lag equal to the difference between the interpolated point and the nearest point.
% simulated_elc.m :  Given a model light curve, generate a list of time-tagged events that their rate follows the model light curve. The events are generated between the first and last time in the model light curve. Note that the first event is always at the first time in the model light curve. The rate at each time is calculated using a second order approximation: Tau0*(1 + dTau/dt + [dTau/dt]^2).
% specwin.m :  Calculate the spectral window of a time series. THIS FUNCTION IS BEING REPLACED BY period.m
% stdfilt1.m :  One dimensional StD filter.
% taper.m :  Generate a taper function for a timeseries. Taper function is a weight/window function in the time domain. This can be used in order to give reduce weight to data found near the edges of a timeseries in order to reduce aliasing with the series length.
