function OutJD=solarlong2jd(SolarLong,ApDate,EquinoxType,AlgoType)
% Time in which the Sun is in a given solar longitude
% Package: celestial.SolarSys
% Description: Calculate time in which the Sun is in a given solar 
%              longitude, using low accuracy formulae (15 min. precision).
% Input  : - Vector of solar longitude (radians).
%          - Matrix of Approximate date of solar longitude [Year, Month].
%          - Equinox for input longitude:
%            'g' - True Geometric referred to the mean equinox of date.
%            'a' - Apparent, referred to the true equinox of date. 
%            'j' - J2000, referred to the J2000.0 equinox. (default).
%          - Algorithm:
%            'l' - low accuracy, default.
%            'e' - Exact solution.
% Output : - Vector of JDs
% See Also : suncoo1.m; suncoo.m
% Reference : Ofek, E. 2000, JIMO, 28, 176
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Sep 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: OutJD=celestial.SolarSys.solarlong2jd(0,[2000 3] ,'a','l')
% Reliable: 2
%--------------------------------------------------------------------
if (nargin==2),
   EquinoxType = 'j';
   AlgoType    = 'l';
elseif (nargin==3),
   AlgoType    = 'l';
elseif (nargin==4),
   % no default
else
   error('Illigal number of input arguments');
end
RAD = 180./pi;
n   = 0.9856474;

Year  = ApDate(:,1);
Month = ApDate(:,2);
Day   = 15.*ones(size(Year));
Frac  = zeros(size(Year));

% start JD
JD    = [celestial.time.julday([Day, Month, Year, Frac])]';


switch AlgoType
 case 'l'
    %--- low accuracy ---
    [RA,Dec,R,SL]=celestial.SolarSys.suncoo(JD,EquinoxType);
    AddDay = n.*(SolarLong - SL);
    
    I = 0;
    while (max(abs(SolarLong-SL))>(0.002./RAD)),
       I = I + 1;
       JD = JD + AddDay;
    
       [RA,Dec,R,SL]=celestial.SolarSys.suncoo(JD,EquinoxType);
       AddDay = n.*(SolarLong - SL);
    end
    OutJD = JD;
 case 'e'
    %--- exact solution ---
    N     = length(SolarLong);
    OutJD = zeros(N,1);
    for I=1:1:N,
       VecJD = [JD(I)-30:0.01:JD(I)+30].';
       [SunRA,SunDec,SunRV,VecSunLong] = celestial.SolarSys.suncoo(VecJD,EquinoxType);
       [Min,MinInd] = min(abs(VecSunLong - SolarLong(I)));
       OutJD(I)     = VecJD(MinInd);
    end

 otherwise
    error('Unknown AlgoType option');
end
