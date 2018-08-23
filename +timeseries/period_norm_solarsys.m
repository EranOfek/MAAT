function PS=period_norm_solarsys(Data,FreqVec,Norm,Delta,R,STO,G,PR)
% Normalized power spectrum for a Solar System object
% Package: timeseries
% Description: Calculate the normalized normal power spectrum of a times
%              series for a solar system object orbiting the Sun. Taking
%              into account ligh travel effect, phase correction due to
%              phase angle, and phase angle changes in brightness.
% Input  : - Two column matrix containing the time series
%            [Time, measurment] or [Time, measurment, error].
%            Where time is measured in days.
%          - Frequency range and interval in which to calculate the
%            power spectrum.
%            This is a column vector of frequencies at which to
%            calculate the power spectrum.
%          - Normalization method:
%            'Var' - Normalize by variance (Default).
%                    If empty use default.
%            'Amp' - Normalize by amplitude.
%          - Observer-Target distance [AU].
%          - Sun-Target distance [AU].
%          - Sun-Target-Observer angle (phase angle) [deg].
%            This should be negative if the target is west of opposition.
%          - Slope parameter in the HG system. Default is 0.15.
%          - Retrograde (-1) or prograde (1) rotation. Default is 1.
% Output : - Two columns matrix of the un-normalized power spectrum
%            [frequency, power_spec].
% See also: period.m
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Jul 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------
RAD = 180./pi;

Def.Norm = 'Var';
Def.G    = 0.15;
Def.PR   = 1;
if (nargin==6)
   G     = Def.G;
   PR    = Def.PR;
elseif (nargin==7)
   PR    = Def.PR;
elseif (nargin==8)
    % do nothing
else
    error('Illegal number of input arguments');
end
if (isempty(Norm))
    Norm = Def.Norm;
end

Col.T = 1;
Col.M = 2;
T       = Data(:,Col.T);
N       = numel(T);
Nf      = numel(FreqVec);

M       = Data(:,Col.M) - mean(Data(:,Col.M));
PS      = zeros(Nf,2);
PS(:,1) = FreqVec;


Phi1 = exp(-3.33.*tan(0.5.*STO./RAD).^0.63);
Phi2 = exp(-1.87.*tan(0.5.*STO./RAD).^1.22);
M    = M - (5.*log10(R.*Delta) - 2.5.*log10((1-G).*Phi1 + G.*Phi2));
M    = M - mean(M);   % subtract mean again

c    = 86400.*get_constant('c')./get_constant('au');  % speed of light [AU/day]

for FreqInd=1:1:Nf
   Tcorr = T - Delta./c;                              % light time correction
   Tcorr = Tcorr + PR.*STO./(360.*FreqVec(FreqInd));  % phase correction
   PS(FreqInd,2) = abs(sum(M.*exp(-2.*pi.*1i.*Tcorr.*FreqVec(FreqInd)))).^2./N;  
end

switch lower(Norm)
 case 'amp'
    % do nothing
 case 'var'
    PS(:,2) = PS(:,2)./var(M);
 otherwise
    error('Unknwon normalization option');
end
