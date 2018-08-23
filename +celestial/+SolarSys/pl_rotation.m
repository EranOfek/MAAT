function [W,Coo,AN]=pl_rotation(Date,Planet)
% Calculate planetary rotation parameters
% Package: celestial.SolarSys
% Description: Calculate planetary rotation parameters.
% Input  : - Date or JD.
%            If one column vector is given then one JD per line.
%            If 3 columns matrix is given then [D M Y] per line. (for UT=0).
%            If 4 columns matrix is given then [D M Y frac_day] per line.
%          - Planet name:
%            'Sun' | 'Mercury' | 'Venus' | 'Earth' | 'Mars' | 'JupiterI' |
%            'JupiterII' | 'JupiterIII' | 'Saturn' |
%            'Uranus' | 'Neptune' | 'Pluto'
% Output : - Location of the prime meridian (W) [radians] measured along
%            the planet's equator in an easterly direction with respect
%            to the planet's north pole from the node (loacated at right
%            asc. pi/2+RA_0) of the planet's equator on the standard
%            equator. If W increases with time, the planet has direct
%            rotation, otherwise its retrograde.
%          - [RA_0, Dec_0], standard equatorial coordinates with equinox
%            J2000 at epoch J2000 in radians.
%          - Aproximate (oscilating)  longitude of ascending node
%            of the planet' orbit referred tomean equinox of date
%            in radians
%            (In the case of the Sun return NaN).
% Notes  : W for Jupiter, Saturn and Uranus refer to the rotation of
%          their magnetic fields (system III). On Jupiter system I refers
%          to the mean atmospheric equatorial rotation. System II referes
%          to the mean atmospheric rotation north of the south component 
%          of the north eq. belt, and south of the north comp. of the
%          south eq. belt.
% Reference: Explanatory supplement P. 705 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [W,Coo,AN]=celestial.SolarSys.pl_rotation([1 1 2000],'JupiterII');
% Reliable: 2
%---------------------------------------------------------------------------

RAD = 180./pi;

SizeDate = size(Date);
Nt       = SizeDate(1);
ColN     = SizeDate(2);

if (ColN==1),
   % already in JD
   JD = Date;
else
   JD = celestial.time.julday(Date).';
end

T = (JD - 2451545.0)./36525;
D = JD - 2451545.0;

switch Planet
 case 'Sun'
    RA_0  = 286.13.*ones(Nt,1);
    Dec_0 = 63.87.*ones(Nt,1);
    W     = 84.10 + 14.1844000.*D;
    AN    = NaN;
 case 'Mercury'
    RA_0  = 281.01 - 0.003.*T;
    Dec_0 = 61.45  - 0.005.*T;
    W     = 329.71 + 6.1385025.*D;
    AN    = 48.330893 + 1.1861890.*T ...
       + 0.00017587.*T.^2 + 0.000000211.*T.^3;
 case 'Venus'
    RA_0  = 272.72.*ones(Nt,1);
    Dec_0 = 67.15.*ones(Nt,1);
    W     = 160.26 - 1.4813596.*D;
    AN    = 76.679920 + 0.9011190.*T ...
       + 0.00040665.*T.^2 - 0.000000080.*T.^3;
 case 'Earth'
    RA_0  = 0.00 - 0.641.*T;
    Dec_0 = 90.00 - 0.557.*T;
    W     = 190.16 + 360.9856235.*D;
    AN    = zeros(size(W));
 case 'Mars'
    RA_0  = 317.681 - 0.108.*T;
    Dec_0 = 52.886 - 0.061.*T;
    W     = 176.868 + 350.8919830.*D;
    AN    = 49.558093 + 0.7720923.*T ...
       + 0.00001605.*T.^2 + 0.000002325.*T.^3;
 case 'JupiterIII'
    RA_0  = 268.05 - 0.009.*T;
    Dec_0 = 64.49 + 0.003.*T;
    W     = 284.95 + 870.5360000.*D;
    AN    = 100.464441 + 1.0209550.*T ...
       + 0.00040117.*T.^2 + 0.000000569.*T.^3;
 case 'JupiterII'
    RA_0  = 268.05 - 0.009.*T;
    Dec_0 = 64.49 + 0.003.*T;
    W     = 43.3 + 870.270.*D;
    AN    = 100.464441 + 1.0209550.*T ...
       + 0.00040117.*T.^2 + 0.000000569.*T.^3;
 case 'JupiterI'
    RA_0  = 268.05 - 0.009.*T;
    Dec_0 = 64.49 + 0.003.*T;
    W     = 67.1 + 877.900.*D;
    AN    = 100.464441 + 1.0209550.*T ...
       + 0.00040117.*T.^2 + 0.000000569.*T.^3;
 case 'Saturn'
    RA_0  = 40.58 - 0.036.*T;
    Dec_0 = 83.54 - 0.004.*T;
    W     = 38.90 + 810.7939024.*D;
    AN    = 113.665524 + 0.8770979.*T ...
       - 0.00012067.*T.^2 - 0.000002380.*T.^3;
 case 'Uranus'
    RA_0  = 257.43.*ones(Nt,1);
    Dec_0 = -15.10.*ones(Nt,1);
    W     = 203.81 - 501.1600928.*D;
    AN    = 74.005947 + 0.5211258.*T ...
       + 0.00133982.*T.^2 + 0.000018516.*T.^3;
 case 'Neptune'
    N     = (359.28 + 54.308.*T)./RAD;
    RA_0  = 299.36 + 0.70.*sin(N);
    Dec_0 = 43.46 - 0.51.*cos(N);
    W     = 253.18 + 536.3128492.*D - 0.48.*sin(N);
    AN    = 131.784057 + 1.1022057.*T ...
       + 0.00026006.*T.^2 - 0.000000636.*T^3;
 case 'Pluto'
    RA_0  = 313.02.*ones(Nt,1);
    Dec_0 = 9.09.*ones(Nt,1);
    W     = 236.77 - 56.3623195.*D;
    AN    = 110.2971389 + (1.396971 - 0.038230).*T;  
 otherwise
    error('Unknown planet');
end

RA_0  = RA_0./RAD;
Dec_0 = Dec_0./RAD;
W     = W./RAD;
AN    = AN./RAD;
Coo   = [RA_0, Dec_0];
