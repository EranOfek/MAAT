function Kcorr=kcorr(Spec,Z1,FiltFam1,FiltName1,MagSys1,Z2,FiltFam2,FiltName2,MagSys2)
% Calculate k-correction.
% Package: AstroUtil.spec
% Description: Calculate k-correction. Given a spectrum two filters, and
%              their redshifts calculate the k-correction of the first
%              filter minus the second filter. This is calculated by
%              shifting the spectrum to z1, measuring the synthetic
%              magnitude in the first filter, then shiting to z2
%              and measuring the synthetic magnitude in the second filter.
% Input  : - Spectrum [Wavelength(Ang), Flux(F_lambda)] of source.
%            Assume the spectrum is at z=0.
%            See convert.flux.m to transform from general units to F_lambda.
%          - Vector of redshift in which to apply the first filter.
%          - First filter family. See get_filter.m for details.
%          - First filter name. See get_filter.m for details.
%          - First magnitude system. See get_filter.m for details.
%          - Vector of redshift in which to apply the second filter.
%            This vector should be of the same length as z1.
%          - Second filter family. See get_filter.m for details.
%          - Second filter name. See get_filter.m for details.
%          - Second magnitude system. See get_filter.m for details.
% Output : - The magnitude difference between the first and second filters,
%            applay to the spectrum at different redshifts.
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Apr 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: Oke & Sandage 1968 ApJ, 154, 21
% Example: Kcorr=AstroUtil.spec.kcorr(Spec,2,'GALEX','NUV','AB',1,'2MASS','J','Vega')
% Reliable: 2
%--------------------------------------------------------------------------

SpecZ1 = AstroUtil.spec.shift_spec(Spec,Z1,'w','f');   % Flux should be corrected by (1+z)
SpecZ2 = AstroUtil.spec.shift_spec(Spec,Z2,'w','f');   % Flux should be corrected by (1+z)

Mag1 = AstroUtil.spec.synphot(SpecZ1,FiltFam1,FiltName1,MagSys1);
Mag2 = AstroUtil.spec.synphot(SpecZ2,FiltFam2,FiltName2,MagSys2);
Kcorr = Mag1 - Mag2;
