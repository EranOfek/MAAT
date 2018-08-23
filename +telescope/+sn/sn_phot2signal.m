function S=sn_phot2signal(SN,Sigma,B,R)
% Convert photometry S/N of a Gaussian source to signal
% Description: Given a target S/N, image background, readnoise and PSF
%              Gaussian sigma, calculate the total count in PSF that
%              will give the target S/N for PSF photometry.
% Input  : - Target S/N.
%          - Sigma of PSF [pixels].
%          - Background [e/pix].
%          - Readout noise [e].
%          - PSF radius [pix]. Default is 20.
% Output : - Source signal [electrons].
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Apr 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=telescope.sn.sn_phot2signal(5,2,300,5)
% Reliable: 2
%--------------------------------------------------------------------------

Thresh = 1e-3;
DeltaS = 0.01;

S1 = 1;
S2 = S1 + DeltaS;
SN1 = telescope.sn.sn_phot_psf(S1,Sigma,B,R);
SN2 = telescope.sn.sn_phot_psf(S2,Sigma,B,R);

while (abs(SN2-SN)>Thresh),
    S1 = S1.*SN./SN1;
    S2 = S1 + DeltaS;
    SN1 = telescope.sn.sn_phot_psf(S1,Sigma,B,R);
    SN2 = telescope.sn.sn_phot_psf(S2,Sigma,B,R);
end
S = S2;
%sn_psf_phot(S,Sigma,B,R)