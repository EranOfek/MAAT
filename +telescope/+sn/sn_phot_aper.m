function SN=sn_phot_aper(S,Sigma,B,R,Radius)
% Calculate the aperture photometry S/N for a Gaussain source
% Package: telescope.sn
% Description: Calculate the S/N (signal-to-noise ratio) for a point
%              source with a symmetric Gaussian profile for aperture
%              photometry.
% Input  : - Source signal [electrons].
%          - Sigma of PSF [pixels].
%          - Background [e/pix].
%          - Readout noise [e].
%          - Aperture photometry radius [pix].
% Output : - Theoretical S/N of aperture photometry.
% See also: sn_psf_phot.m, optimal_phot_aperture.m
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Apr 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: SN=telescope.sn.sn_phot_aper(1000,2,100,10,5)
% Reliable: 2
%--------------------------------------------------------------------------

Phi = 1 - exp(-Radius.^2./(2.*Sigma.^2));
SN = S.*Phi./sqrt(S.*Phi + pi.*Radius.^2.*(B+R.^2));
