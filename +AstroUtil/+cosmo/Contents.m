%
% Contents file for package: cosmology
% Created: 29-Dec-2015
%---------
% ad_dist.m :  Calculate the filled beam angular diameter distance between two redshifts along the line of sight.
% ad_q_dist.m :  Compute filled-beam angular-diameter distance to an object in a flat Universe with constant equation of state (p=w\rho; i.e., quintessence).
% cdt_dz.m :  Calculate the differential cdt/dz in the FLRW geometry.
% comoving_dist.m :  Calculate the line of sight comoving distance.
% comoving_volume.m :  Calculate the differential comoving volume (d\Omega dz) at redshift z and the total comoving volume from redshift 0 to z.
% cosmo_pars.m :  Return the cosmological parameters as measured by various experiments.
% crit_surface_density.m :  Calculates the critical surface density for gravitational lensing, given the cosmology and redshifts of the lens and source.
% delta_vir_z.m :  Calculate the virial overdensity \Delta_{vir}, as a function of redshift z, and cosmological paramaeters.
% dist_mod2dist.m :  Convert distance modulous to luminosity distance and redshift.
% e_z.m :  Calculate E(z) cosmological function, which is proportional to the time derivative of the logarithm of the scale factor.
% hubble_z.m :  Compute the Hubble parameter as a function of redshift. (Assuming matter dominated universe - Z<1000).
% inv_comoving_volume.m :  Use the cosmological volume to calculate the corresponding redshift.
% inv_e_z.m :  Calculate 1/E(z) cosmological function, in which E(z) is proportional to the time derivative of the logarithm of the scale factor.
% inv_lum_dist.m :  Given the distance modulus, use the luminosity distance to calculate the corresponding redshift.
% lookback_time.m :  Compute the cosmological lookback time, between two events in redshift z1 and z2, and given the cosmology. (Assuming matter dominated universe - Z<1000).
% lum_dist.m :  Compute luminosity distance from redshift and cosmological parameters. Given the object spectra, calculate also the K-correction.
% matter_density.m :  Calculate the mean matter density in the Universe.
% omega_m_lambda_lines.m :  Given a universe with \Omega_{m} and \Omega_{\Lambda} contributions, and given \Omega_{m} vector, find for each value of \Omega_{m}: (i) the value of \Omega_{\Lambda} for which the universe will expand forever; (ii) The \Omega_{\Lambda} criterion for which there have been no singularity in the past (rather than Big Bang its early history consisted of a period of gradually slowing contraction to a minimum radius before begining its current expansion).
% omega_z.m :  Calculate \Omega_{m} as a function of redshift z.
% tran_comoving_dist.m :  Calculate the transverse comoving distance. Given the transverse comoving distance (D_M), the comoving distance between two events at the same redshift, but separated on the sky by some angle \delta\theta is: D_M\delta\theta.
