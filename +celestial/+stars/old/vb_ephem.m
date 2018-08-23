function [D,PA,RV]=vb_ephem(Time,OrbElem)
%------------------------------------------------------------------------------
% vb_ephem function                                                      ephem
% Description: Given orbital elements of a visual binary, calculate its
%              ephemeris in a give dates.
% Input  : - Vector of JD.
%          - Vector of visual binary orbital elements:
%            [T,   a,     e, Om, w,  i,  P], units are:
%            [JD,  arcsec,[],deg,deg,deg,Julian_yr].
%            Alternatively, this can be a structure, containing the
%            fields .T, .a, .e, .Om, .w, .i, .P.
% Output : - Vector of distances between binary star [arcsec], per
%            each time.
%          - Vector of position angle [radians], per each time. Referred,
%            to the Equinox in which "Om" (long. of ascending node) is
%            given.
%          - Vector of [parallax * radial_velocity], in [arcsec km/s],
%            assuming negligible secondary mass.
% Tested : Matlab 4.2
%     By : Eran O. Ofek                  November 1995
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Web example: http://astroclub.tau.ac.il/ephem/VisualDoubleStars/
% Reliable: 2
%------------------------------------------------------------------------------
RAD        = 180./pi;
JulianYear = 365.25;
DAY_SEC    = 86400;
AU_KM      = 149597870;    % km

if (isstruct(OrbElem)==1),
   % do nothing - use structure as is
   El = OrbElem;
else
   % set structure
   El.T = OrbElem(1);
   El.a = OrbElem(2);
   El.e = OrbElem(3);
   El.w = OrbElem(4);
   El.i = OrbElem(5);
   El.P = OrbElem(6);
end
% convert period from years to days (assuming period is in Julian years).
El.P = El.P.*JulianYear;

% convert angles to radians
El.w   = El.w./RAD;
El.i   = El.i./RAD;
El.Om  = El.Om./RAD;

El.q   = El.a.*(1 - El.e);
N      = 2.*pi./El.P;    % mean motion [rad/day]
% call kepler_elliptic with the mean anomaly option:
[Nu,R] = kepler_elliptic(N.*(Time-El.T), El.q, El.e, NaN);   

PA = El.Om + atan2( sin(Nu + El.w).*cos(El.i), cos(Nu + El.w) );
PA = (PA./(2.*pi) - floor(PA./(2.*pi))).*2.*pi;   % range [0-2*pi]
D  = R.*cos(Nu + El.w)./cos(PA - El.Om);


K2 = N.*El.a.*sin(El.i)./sqrt(1 - El.e.^2);   % Secon. RV amplitude
                                              % [parallax*au/day]
RV = K2.*(El.e.*cos(El.w) + cos(Nu + El.w));  % Secon. RV 

% convert to [parallax * km/s]
RV = RV.* AU_KM./DAY_SEC;




