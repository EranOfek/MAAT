function PeakWave=wein(Temp,UnitsTemp)
% Wein law
% Package: AstroUtil.spec
% Description: Apply Wein law - return the peak wavelength of a
%              black body at a given temperature.
% Input  : - Temperature.
%          - Units of temperature, 
%            {'erg'|'J'|'Hz'|'A'|'eV'|'T'|'me'|'mp'|'cal'|'Btu'|
%             'kWh'|'TNT','gr'}, default is 'T' (Kelvin),
%            see convert_energy.m for options.
% Output : - Wavelength of peak intensity of black body spectrum [Ang].
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Feb 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: PeakWave=wein(1,'eV'); PeakWave=wein(5770,'T');
% Reliable: 1
%--------------------------------------------------------------------------

DefUnitsTemp = 'T';
if (nargin==1)
   UnitsTemp = DefUnitsTemp;
elseif (nargin==2)
   % do nothing
else
   error('Illegal number of input arguments');
end

TempK    = convert.energy(UnitsTemp,'T',Temp);
PeakWave = 2.898e7./TempK;   % [Ang]

