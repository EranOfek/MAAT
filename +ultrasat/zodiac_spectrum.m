function OutSpec=zodiac_spectrum(varargin)
% Get the Zodiac light spectrum
% Package: @AstSpec
% Description: Return the zodiac spectrum as adopted from the HST STIS
%              handbook. The high zodiacal ligh is defined where V=22.1
%              mag/arcsec^-2.
% Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Wave' - Vector of wavelength [Ang] in which to calculate
%                   spectrum. If empty, then use original data.
%                   Default is empty.
%            'BackType' - Background type. Default is 'zodi'.
%            'OutType'  - 'mat'|'astspec'. Default is 'mat'.
%            'InterpMethod' - Default is 'linear'.
% Output : - Zodiacal ligh spectrum
%            [wavelength(Ang), Flux(erg/cm^2/s/A/arcsec^2)]
% Reference: https://hst-docs.stsci.edu/display/STISIHB/6.6+Tabular+Sky+Backgrounds
%            but there is a discrepency with:
%            http://www.stsci.edu/hst/wfc3/design/documents/handbooks/currentIHB/c09_exposuretime08.html#389841
%            According to the HST help desk the STIS table should be used.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Nov 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=ultrasat.zodiac_spectrum;
%          % to verify normalization: synphot(Spec,'Johnson','V','Vega')
%          S=ultrasat.zodiac_spectrum;
% Reliable: 1
%--------------------------------------------------------------------------


DefV.Wave                 = [];
DefV.BackType             = 'zodi';   % 'zodi' | 'earthshine' | 'total' | 'all'
DefV.OutType              = 'mat';
DefV.InterpMethod         = 'linear'; 
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

ColWave = 1;
switch lower(InPar.BackType)
    case 'zodi'
        Col = 3;
    case 'earthshine'
        Col = 2;
    case 'total'
        Col = 4;
    case 'all'
        Col = [2 3 4];
    otherwise
        error('Unknown BackType option');
end

InterpMethod = 'linear';



% Wave(Ang), High Earthshine, High zodi, Total back
% erg/s/cm^2/A/arcsec^2
Spec=[
1000 2.41e-23 9.69e-29 2.41e-23
1100 4.38e-22 1.04e-26 4.38e-22
1200 4.01e-23 1.08e-25 4.03e-23
1300 7.41e-25 6.59e-25 1.40e-24
1400 4.29e-25 2.55e-24 2.98e-24
1500 4.16e-25 9.73e-24 1.01e-23
1600 2.55e-25 2.35e-22 2.35e-22
1700 7.89e-25 7.21e-21 7.21e-21
1800 9.33e-23 1.53e-20 1.54e-20
1900 4.39e-22 2.25e-20 2.29e-20
2000 1.01e-21 3.58e-20 3.68e-20
2100 1.60e-21 1.23e-19 1.25e-19
2200 7.49e-22 2.21e-19 2.22e-19
2300 3.32e-22 1.81e-19 1.81e-19
2400 2.50e-22 1.83e-19 1.83e-19
2500 2.39e-22 2.53e-19 2.53e-19
2600 5.62e-22 3.06e-19 3.06e-19
2700 6.77e-21 1.01e-18 1.02e-18
2800 2.03e-21 2.88e-19 2.90e-19
2900 4.32e-20 2.08e-18 2.12e-18
3000 9.34e-20 1.25e-18 1.35e-18
3100 2.07e-19 1.50e-18 1.70e-18
3200 3.60e-19 2.30e-18 2.66e-18
3300 4.27e-19 2.95e-18 3.38e-18
3400 6.40e-19 2.86e-18 3.50e-18
3500 8.20e-19 2.79e-18 3.61e-18
3600 1.06e-18 2.74e-18 3.80e-18
3700 1.22e-18 3.32e-18 4.54e-18
3800 1.23e-18 3.12e-18 4.35e-18
3900 1.52e-18 3.34e-18 4.86e-18
4000 2.38e-18 4.64e-18 7.01e-18
4250 2.38e-18 4.65e-18 7.03e-18
4500 2.86e-18 5.58e-18 8.44e-18
4750 2.79e-18 5.46e-18 8.25e-18
5000 2.63e-18 5.15e-18 7.77e-18
5250 2.67e-18 5.37e-18 8.04e-18
5500 2.58e-18 5.34e-18 7.92e-18
5750 2.54e-18 5.40e-18 7.94e-18
6000 2.42e-18 5.25e-18 7.67e-18
6250 2.26e-18 5.02e-18 7.28e-18
6500 2.17e-18 4.92e-18 7.09e-18
6750 2.07e-18 4.79e-18 6.87e-18
7000 1.93e-18 4.55e-18 6.48e-18
7250 1.85e-18 4.43e-18 6.29e-18
7500 1.74e-18 4.23e-18 5.97e-18
7750 1.63e-18 4.04e-18 5.67e-18
8000 1.56e-18 3.92e-18 5.49e-18
8250 1.48e-18 3.76e-18 5.23e-18
8500 1.35e-18 3.50e-18 4.85e-18
8750 1.31e-18 3.43e-18 4.74e-18
9000 1.22e-18 3.23e-18 4.44e-18
9250 1.15e-18 3.07e-18 4.21e-18
9500 1.10e-18 2.98e-18 4.08e-18
9750 1.04e-18 2.86e-18 3.91e-18
10000 1.00e-18 2.78e-18 3.78e-18
10250 9.54e-19 2.67e-18 3.63e-18
10500 9.04e-19 2.56e-18 3.46e-18
10750 8.41e-19 2.41e-18 3.25e-18
11000 8.03e-19 2.31e-18 3.11e-18];

% selected columns
SpecInt = Spec(:,Col);


if isempty(InPar.Wave)
    % use default wavelength
    Wave = Spec(:,ColWave);
else
    Wave = InPar.Wave;
    SpecInt = interp1(Spec(:,ColWave),SpecInt,Wave,InPar.InterpMethod);
end

switch lower(InPar.OutType)
    case 'mat'
        OutSpec = [Wave, SpecInt];
    case 'astspec'
        OutSpec = AstSpec;
        OutSpec.Wave = Wave;
        OutSpec.Int  = SpecInt;
        OutSpec.WaveUnits = 'Ang';
        OutSpec.IntUnits  = 'erg*cm^-2*s^-1*Ang^-1*arcsec^-2';
        OutSpec.ObjName   = 'Zodiac spectrum';
        OutSpec.source    = 'HST STIS handbook';
        OutSpec.z         = 0;
    otherwise
        error('Unknwon OutType option');
end


