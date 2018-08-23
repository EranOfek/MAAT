%
% Contents file for package: AstroStat
% Created: 29-Dec-2015
%---------
% bc_a.m :  Bias correction and acceleration for bootstrap and Jackknife estimetors of the confidence interval. The bc_a bias coorection and acceleration can be used to estimate the bias corrected and accelerated confiedence interval (CI).
% bootstrap_std.m :  Given an estimator (given by a function) calculate the Bootstrap StD for this estimator.
% cel_coo_rnd.m :  Generate random coordinates on the celestial sphere. The program can applay matrix rotation, reject coordinates, and generate non uniform coordinates. [preliminary version].
% confint_probdist.m :  Calculate two-sided confidence interval of a given numerical probability distribution function.
% corrsim.m :  Calculate the correlation between two vectors and use the bootstrap method to estimate the probability to get a correlation larger than the observed correlation. The function ignores NaN values.
% corrsim_cov.m :  Given a matrix with N columns, calculate the correlation between each pair of columns and use the bootstrap method to estimate the probability to get a correlation larger than the observed correlation. The function ignores NaN values.
% err_cl.m :  Given a vector of data points, calculate the lower and upper bounds of an interval that contains a given precentage (P) of the data. (1-P)/2 of the data points are below and above the lower and upper bounds, respectively. Remove NaNs prior to calculation.
% error2ensemble.m :  Generate a realization of data points given the data probability distribution function.
% hist2d.m :  calculate the 2-D histogram of 2-D data set.
% jackknife.m :  Given an estimator (given by a function), calculate the Jackknife StD and the first order Quenouille-Tukey Jacknife bias for this estimator. Notes: - the estimator should be continues. - if Bias/StD~<0.25 then the bias is probably not an issue. - The bias estimate is not reliable if the estimator is an unsmooth statistic (e.g., median).
% max_likelihood.m :  Given a numerical probability distribution and list of 'events', calculate the 'likelihood' for the events given the probability distribution. In addition, calculate the ML ratio test and preform Monte-Carlo simulations by generating realizations of events given the probability distribution and calculate the likelihood-probability distribution.
% maxnd.m :  Return the global maximum of a N-D matrix and its indices. This is equivalent to max(Data(:)), but it also returns the indices of the global maximum.
% mean_error.m :  Calculate the error on the mean using std/sqrt(N).
% meannd.m :  Return the global mean of a N-D matrix.
% mediannd.m :  Return the global median of a N-D matrix.
% minnd.m :  Return the global minimum of a N-D matrix and its indices. This is equivalent to min(Data(:)), but it also returns the indices of the global minimum.
% mode_fit.m :  Estimate the mode of an array by fitting a Gaussian to the histogram of the array around its median. Return also the Sigma of the Gaussian fit.
% moment_2d.m :  Calculate first and second moments of a 2D matrix.
% nanrstd.m :  Robust nanstd. Estimating the std (like nanstd.m), based on the 50-percentile of the distribution.
% noiser.m :  Add multiple noise components to a vector.
% poissconf.m :  Given the number of observed events, calculates the two sided upper and lower confidence intervals, assuming Poisson statistics. Below N=140 use the Gehrels (1986) algorithm. Above this, use sqrt(N) approximation.
% prob2find_inr.m :  Given a density (number per unit area), and a distance from a point, calculate the the probability to find an object within radius R from the point (assuming Poisson statistics).
% psigma.m :  Return the two sided or one sided probability for a given sigma level.
% rand_circle.m :  Generate random number equally distributed inside a unit circle.
% rand_ps.m :  Generate a random realization of a time series with a given power spectrum (e.g., power-law) and optional gaussian measuments errors.
% rand_range.m :  Generate uniformly random number in a given range. The numbers can be uniform either in linear space or log10 space.
% randgen.m :  Random numbers generator for arbitrary function given by a numerical distribution.
% randinpolygon.m :  Generate random positions inside a polygon, defined on a plane or a sphere.
% rangend.m :  Return the global Range of a N-D matrix.
% realhist.m :  Calculate histogram for a dataset in a given range.
% stdnd.m :  Return the global StD of a N-D matrix.
% sumnd.m :  Return the global sum of a N-D matrix.
% symerror.m :  Given a symbolic expression and the variables in the expression, calculate the symbolic error function of the expression, with respect to the variables. The output error expression contains an error variable named D_"original_var" for each variable in the input expression.
% symerror_calc.m :  Given a symbolic expression, the names of the variables in the expression and the value of the variables and errors, calculate the symbolic error function and evaluate it.
% wmean.m :  Calculated the weighted mean of a sample (ignoring NaNs).
% xcorr2p.m :  Given two lists of coordinates (selected in the same region), calculate the cross correlation function as a function of distance. The cross-correlation can be calculated using the DR, RR or DR2 methods.
