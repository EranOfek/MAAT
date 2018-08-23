function S=sn_det2signal(SN,Sigma,B,R)
% Convert detection S/N of a Gaussian PSF to signal
% Package: telescope.sn
% Description: Given detection S/N calculate the PSF signal.
% Input  : - S/N for detection
%          - Width of Gaussian PSF (sigma in pix).
%          - Background [e-]
%          - Readnoise [e-].
% Output : - Signal [e-]
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Sep 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=telescope.sn.sn_det2signal(10,2,100,0)
% Reliable: 2
%--------------------------------------------------------------------------

S  = SN.*sqrt(4.*pi.*(B+R.^2).*Sigma.^2);
