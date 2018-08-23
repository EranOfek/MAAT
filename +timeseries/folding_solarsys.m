function [Phase,Mag]=folding_solarsys(Data,Period,Delta,R,STO,G,PR);
% Folding a timeseries for a solar system object
% Package: timeseries
% Description: Folding a time series of a solar system object into a
%              period, and calculate the phase for each data point in
%              the time series.  Taking
%              into account ligh travel effect, phase correction due to
%              phase angle, and phase angle changes in brightness.
% Input  : - Matrix of [Time, Mag, Err], where time is measured in days.
%          - Period to fold into [day].
%          - Observer-Target distance [AU].
%          - Sun-Target distance [AU].
%          - Sun-Target-Observer angle (phase angle) [deg].
%            This should be negative if the target is west of opposition.
%          - Slope parameter in the HG system. Default is 0.15.
%          - Retrograde (-1) or prograde (1) rotation. Default is 1.
% Output : - Vector of phases (0 to 1) corresponding to each line in the
%            first input argument.
%          - Vector of minor planet corrected HG magnitudes.
% Tested : Matlab 3.5
%     By : Eran O. Ofek                    Jul 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------
RAD = 180./pi;

Def.G    = 0.15;
Def.PR   = 1;
if (nargin==5),
   G     = Def.G;
   PR    = Def.PR;
elseif (nargin==6),
   PR    = Def.PR;
elseif (nargin==7),
    % do nothing
else
    error('Illegal number of input arguments');
end

Time = Data(:,1);
Mag  = Data(:,2);
Err  = Data(:,3);

Phi1 = exp(-3.33.*tan(0.5.*STO./RAD).^0.63);
Phi2 = exp(-1.87.*tan(0.5.*STO./RAD).^1.22);
Mag  = Mag - (5.*log10(R.*Delta) - 2.5.*log10((1-G).*Phi1 + G.*Phi2));

c    = 86400.*constant.c./get_constant.au;  % speed of light [AU/day]

Tcorr = Time  - Delta./c;                          % light time correction
Tcorr = Tcorr + PR.*STO.*Period./360;  % phase correction
Phase = Tcorr./Period - floor(Tcorr./Period);
