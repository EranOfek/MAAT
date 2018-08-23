function EarthVel=earth_vel_ron_vondrak(JD,Units)
% Earth barycentric velocity
% Package: celestial.SolarSys
% Description: Calculate the Earth barycentric velocity in respect to the
%              mean equator and equinox of J2000.0, using a version
%              of the Ron & Vondrak (1986) trigonometric series.
% Input  : - Vector of JD.
%          - Output units. Options are:
%            'cgs'  - for [cm].
%            'SI'   - for [m].
%            'agd'  - for [AU] - default.
% Output : - Earth barycentric velocity vector refereed to the equatorial
%            system and true equinox of date. This is a three column
%            matrix [V_X, V_Y, V_Z] in which each row corresponds to one
%            epoch.
% Tested : matlab 7.10
%     By : Eran O. Ofek                    Oct 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: EarthVel=celestial.SolarSys.earth_vel_ron_vondrak(2451545)
% Reliable: 2
%--------------------------------------------------------------------------

Def.Units = 'agd';
if (nargin==1),
  Units = Def.Units;
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (size(JD,2)==1),
   % convert JD to row vector
   JD = JD.';
end

T = (JD - 2451545.0)./36525.0;

% all expressions are given in radians:
%L1  = 4.4026088 + 2608.7903142.*T;
L2  = 3.1761467 + 1021.3285546.*T;
L3  = 1.7534703 +  628.3075849.*T;
L4  = 6.2034809 +  334.0612431.*T;
L5  = 0.5995465 +   52.9690965.*T;
L6  = 0.8740168 +   21.3299095.*T;
L7  = 5.4812939 +    7.4781599.*T;
L8  = 5.3118863 +    3.8133036.*T;
%W   = 3.8103444 + 8399.6847337.*T;
D   = 5.1984667 + 7771.3771486.*T;
Mt  = 2.3555559 + 8328.6914289.*T;
Lt  = 3.8103444 + 8399.6847337.*T;
F   = 1.6279052 + 8433.4661601.*T;


Arg = [L3
       2.*L3 
       L5
       Lt
       3.*L3
       L6
       F
       Lt+Mt
       2.*L5
       2.*L3-L5
       3.*L3-8.*L4+3.*L5
       5.*L3-8.*L4+3.*L5
       2.*L2-L3
       L2
       L7
       L3-2.*L5
       L8
       L3+L5
       2.*L2-2.*L3
       L3-L5
       4.*L3
       3.*L3-2.*L5
       L2-2.*L3
       2.*L2-3.*L3
       2.*L6
       2.*L2-4.*L3
       3.*L3-2.*L4
       Lt+2.*D-Mt
       8.*L2-12.*L3
       8.*L2-14.*L3
       2.*L4
       3.*L2-4.*L3
       2.*L3-2.*L5
       3.*L2-3.*L3
       2.*L3-2.*L4
       Lt-2.*D];

ONES  = ones(size(T));
SinX = [-1719914-2.*T
	 6434+141.*T
	 715  .* ONES
	 715  .* ONES
         486-5.*T
         159  .*ONES
	 0  .*ONES
         39  .*ONES
         33  .*ONES
         31  .*ONES
         8  .*ONES
         8  .*ONES
         21  .*ONES
        -19  .*ONES
         17  .*ONES
         16  .*ONES
         16  .*ONES
         11  .*ONES
         0  .*ONES
        -11  .*ONES
        -7  .*ONES
        -10  .*ONES
        -9  .*ONES
        -9  .*ONES
         0  .*ONES
         0  .*ONES
         8  .*ONES
         8  .*ONES
        -4  .*ONES
        -4  .*ONES
        -6  .*ONES
        -1  .*ONES
	 4  .*ONES
         0  .*ONES
         5  .*ONES
	 5  .*ONES];

CosX = [-25.*ONES
        28007-107.*T
        0 .*ONES
        0 .*ONES
       -236-4.*T
        0 .*ONES
        0 .*ONES
        0 .*ONES
        -10 .*ONES
        1 .*ONES
        -28 .*ONES
        -28 .*ONES
        0 .*ONES
        0 .*ONES
        0 .*ONES
        0 .*ONES
        0 .*ONES
        -1 .*ONES
        -11 .*ONES
        -2 .*ONES
        -8 .*ONES
        0 .*ONES
        0 .*ONES
        0 .*ONES
        -9 .*ONES
        -9 .*ONES
        0 .*ONES
        0 .*ONES
        -7 .*ONES
        -7 .*ONES
        -5 .*ONES
        -1 .*ONES
        -6 .*ONES
        -7 .*ONES
        -5 .*ONES
         0.*ONES];




SinY = [25-13.*T
        25697-95.*T
        6  .*ONES
        0  .*ONES
       -216-4.*T
        2  .*ONES
        0  .*ONES
        0  .*ONES
        -9  .*ONES
        1  .*ONES
        25  .*ONES
        -25  .*ONES
        0  .*ONES
        0  .*ONES
        0  .*ONES
        0  .*ONES
        1  .*ONES
        -1  .*ONES
        -10  .*ONES
        -2  .*ONES
        -8  .*ONES
        0  .*ONES
        0  .*ONES
        0  .*ONES
        -8  .*ONES
        8  .*ONES
        0  .*ONES
        0  .*ONES
        -6  .*ONES
        6  .*ONES
        -4  .*ONES
        -2  .*ONES
        -5  .*ONES
        -6  .*ONES
        -4  .*ONES
	 0  .*ONES];


CosY = [1578089+156.*T
        -5904-130.*T
        -657  .*ONES
        -656  .*ONES
        -446+5.*T
        -147  .*ONES
        26  .*ONES
        -36  .*ONES
        -30  .*ONES
        -28  .*ONES
        8  .*ONES
        -8  .*ONES
        -19  .*ONES
        17  .*ONES
        -16  .*ONES
        15  .*ONES
        -15  .*ONES
        -10  .*ONES
        0  .*ONES
        9  .*ONES
        6  .*ONES
        9  .*ONES
        -9  .*ONES
        -8  .*ONES
        0  .*ONES
        0  .*ONES
        -8  .*ONES
        -7  .*ONES
        4  .*ONES
        -4  .*ONES
        5  .*ONES
        -7  .*ONES
        -4  .*ONES
         0  .*ONES
        -5  .*ONES
        -5  .*ONES];


SinZ = [10+32.*T
        11141-48.*T
        -15  .*ONES
        0  .*ONES
        -94  .*ONES
        -6  .*ONES
        0  .*ONES
        0  .*ONES
        -5  .*ONES
        0  .*ONES
        11  .*ONES
        -11  .*ONES
        0  .*ONES
        0  .*ONES
        0  .*ONES
        1  .*ONES
        -3  .*ONES
        -1  .*ONES
        -4  .*ONES
        -1  .*ONES
        -3  .*ONES
        0  .*ONES
        0  .*ONES
        0  .*ONES
        -3  .*ONES
        3  .*ONES
        0  .*ONES
        0  .*ONES
        -3  .*ONES
        3  .*ONES
        -2  .*ONES
        1  .*ONES
        -2  .*ONES
        -3  .*ONES
        -2  .*ONES
	 0  .*ONES];


CosZ = [ 684185-358.*T
        -2559-55.*T
        -282  .*ONES
        -285  .*ONES
       -193  .*ONES
        -61  .*ONES
        -59  .*ONES
        -16  .*ONES
        -13  .*ONES
        -12  .*ONES
        3  .*ONES
        -3  .*ONES
        -8  .*ONES
        8  .*ONES
        -7  .*ONES
        7  .*ONES
        -6  .*ONES
        -5  .*ONES
        0  .*ONES
        4  .*ONES
        3  .*ONES
        4  .*ONES
        -4  .*ONES
        -4  .*ONES
        0  .*ONES
        0  .*ONES
        -3  .*ONES
        -3  .*ONES
        2  .*ONES
        -2  .*ONES
        2  .*ONES
        -4  .*ONES
        -2  .*ONES
        0  .*ONES
        -2  .*ONES
	-2  .*ONES];

V_X = (sum(SinX.*sin(Arg)) + sum(CosX.*cos(Arg))).';
V_Y = (sum(SinY.*sin(Arg)) + sum(CosY.*cos(Arg))).';
V_Z = (sum(SinZ.*sin(Arg)) + sum(CosZ.*cos(Arg))).';

EarthVel = [V_X, V_Y, V_Z].*1e-8;  % [au/day]

switch lower(Units)
 case 'agd'
    % do nothing - already in [au/day]
 case 'si'
    EarthVel = EarthVel.*convert.units('au/day','m/s');
 case 'cgs'
    EarthVel = EarthVel.*convert.units('au/day','cm/s');
 otherwise
    error('Unknown Units option');
end

