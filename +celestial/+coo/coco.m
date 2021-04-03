function [OutList,TotRot]=coco(InList,InCooType,OutCooType,InUnits,OutUnits)
% Convert between different coordinates (OBSOLETE: use convert_coo)
% Package: celestial.coo
% Description: General coordinate convertor. Convert/precess
%              coordinate from/to Equatorial/galactic/Ecliptic
%              coordinates.
% Input  : - Matrix of input coordinates.
%            First column for Long/RA and second column for lat/Dec.
%          - Type of input coordinates.
%            'j####.#' - equatorial, mean of date
%            (default 'j2000.0'). [Julian].
%            'J####.#' - equatorial, true of date
%            (default 'j2000.0'). [Julian].
%            'g' - galactic.
%            'S' - Super galactic.
%            'c' - CMB dipole (see rotm_coo.m for details)
%            'e' - Ecliptic with J2000.0 equinox.
%          - Type of outpt coordinates.
%            'j####.#' - equatorial, mean of date
%            (default 'j2000.0'). [Julian].
%            'J####.#' - equatorial, true of date
%            (default 'j2000.0'). [Julian].
%            'g' - galactic. (default)
%            'S' - Super galactic.
%            'c' - CMB dipole (see rotm_coo.m for details)
%            'e' - Ecliptic with J2000.0 equinox.
%          - Units for input coordinates.
%            'r' - radians. (default)
%            'd' - degrees.
%            'h' - hours/deg.
%          - Units for outpu coordinates.
%            'r' - radians. (default)
%            'd' - degrees.
%            'h' - hours/deg.
% Output : - Matrix of output coordinates.
%          - Total rotation matrix.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Feb 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: C=celestial.coo.coco(rand(10,2),'j2000.0','e')
% Reliable: 1
%--------------------------------------------------------------------------
import celestial.coo.*

RADIAN = 180./pi;

if (nargin==1)
   InCooType  = 'j2000.0';
   OutCooType = 'g';
   InUnits    = 'r';
   OutUnits   = 'r';
elseif (nargin==2)
   OutCooType = 'g';
   InUnits    = 'r';
   OutUnits   = 'r';
elseif (nargin==3)
   InUnits    = 'r';
   OutUnits   = 'r';
elseif (nargin==4)
   OutUnits   = 'r';
elseif (nargin==5)
   % do nothing
else
   error('Illigal number of input arguments');
end


LenInType  = length(InCooType);
LenOutType = length(OutCooType);

if (ischar(InCooType)==1)
   if (LenInType>1)
      InEquinox  = str2double(InCooType(2:LenInType));
      InEquinoxJD = 2451545.5 + 365.25.*(InEquinox - 2000);
      InCooType  = InCooType(1);
   end
else
   InEquinox = 0;
   InEquinoxJD = InCooType;
   InCooType = 'j';
   InCooType  = InCooType(1);
end

if (ischar(OutCooType)==1)
   if (LenOutType>1)
      OutEquinox = str2num(OutCooType(2:LenOutType));
      OutEquinoxJD = 2451545.5 + 365.25.*(OutEquinox - 2000);
      OutCooType  = OutCooType(1);
   end
else
   OutEquinox = 0;
   OutEquinoxJD = OutCooType;
   OutCooType = 'j';
   OutCooType  = OutCooType(1);
end

InCoo   = zeros(size(InList));
OutCoo  = zeros(size(InList));
OutList = zeros(size(InList));

switch InUnits
 case {'r'}
    InCoo = InList;
 case {'d'}
    % convert deg. to radians
    InCoo = InList./RADIAN;
 case {'h'}
    % convert h/d to radians
    InCoo(:,1) = InList(:,1).*15./RADIAN;
    InCoo(:,2) = InList(:,2)./RADIAN;
 otherwise
    error('Unknown type of input units');
end

% convert coordinates to direction cosines
InCosDir = cosined(InCoo);


RotM1 = diag([1 1 1]);

% calculate the first rotation matrix
switch InCooType
 case {'j','J'}
    if (InEquinox~=2000.0)
       % precess coordinates to J2000.0
       switch InCooType
       case {'j'}
           % mean equinox ...
           RotM1 = rotm_coo('p',InEquinoxJD);
        case {'J'}
           % true equinox ...
           RotM1 = rotm_coo('pd',InEquinoxJD);
        otherwise
           error('Illegal InCooType');
       end
    end
 case {'g'}
    % convert to Equatorial J2000.0
    RotM1 = rotm_coo('G',2451545.5);
 case {'S'}
    % convert to Equatorial J2000.0
    RotM1 = rotm_coo('G',2451545.5)*rotm_coo('SGg');
 case {'c'}
    % convert to Equatorial
    error('not implemented')
 case {'e'}
    % convert to Equatorial J2000.0
    RotM1 = rotm_coo('E',2451545.5);
 otherwise
    error('Unknown input coordinaytes type');
end

RotM2 = diag([1 1 1]);
% calculate the second rotation matrix
switch OutCooType
 case {'j','J'}
    if (OutEquinox~=2000.0)
       % precess coordinates from J2000.0
       switch OutCooType
        case {'j'}   
           % mean equinox ...
           RotM2 = rotm_coo('P',OutEquinoxJD);
        case {'J'}   
           % true equinox ...
           RotM2 = rotm_coo('Pd',OutEquinoxJD);
        otherwise
           error('Illegal OutCooType');
       end           
    end
 case {'g'}
    % convert to galactic
    RotM2 = rotm_coo('g',2451545.5);
 case {'S'}
    % convert to Super galactic
    RotM2 = rotm_coo('gSG')*rotm_coo('g',2451545.5);
 case {'c'}
    % convert to CMB dipole
    error('not implemented')
 case {'e'}
    % convert to ecliptic
    RotM2 = rotm_coo('e',2451545.5);
 otherwise
    error('Unknown output coordinaytes type');
end

% rotate coordinates
TotRot = RotM2*RotM1;
OutCosDir = TotRot*[InCosDir.'];


% convert coordinates from direction cosines
OutCoo = cosined([OutCosDir.']);


I0 = find(OutCoo(:,1)<0);
OutCoo(I0,1) = 2.*pi + OutCoo(I0,1);

switch OutUnits
 case {'r'}
    OutList = OutCoo;
 case {'d'}
    % convert radians to deg.
    OutList = OutCoo.*RADIAN;
 case {'h'}
    % convert radians to h/d
    OutList(:,1) = OutCoo(:,1).*RADIAN./15;
    OutList(:,2) = OutCoo(:,2).*RADIAN;
 otherwise
    error('Unknown type of output units');
end





