function SN=sn_phot_psfn(Flux,PSF,Variance)
% Calculate photometry S/N of a numerical PSF
% Description: Calculate the S/N (signal-to-noise ratio) for a numerical
%              PSF (optimal) photometry.
% Input  : - Flux normalization.
%          - PSF stamp normalized to 1.
%            Alternatively, this can be a cell array of PSFs.
%          - The background variance (including background and readnoise).
%            If PSF is a cell array of PSFs, this need to be a cell array
%            of variances of the same length.
% Output : - Theoretical S/N of PSF photometry.
% See also: sn_psf_phot.m
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Apr 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: SN=telescope.sn.sn_phot_psfn(1000,VecPSF{1},VecStd(1).^2)
% Reliable: 2
%--------------------------------------------------------------------------

if (isnumeric(PSF)),
    PSF = {PSF};
end
if (isnumeric(Variance)),
    Variance = {Variance};
end
Npsf = numel(PSF);
SN   = zeros(Npsf,1);
for Ipsf=1:1:Npsf,
    SN(Ipsf) = sqrt(sum((Flux(Ipsf).*PSF{Ipsf}(:)).^2./(Flux(Ipsf).*PSF{Ipsf}(:) + Variance{Ipsf}(:))));   % sqrt(sum(E^2/V))
end
