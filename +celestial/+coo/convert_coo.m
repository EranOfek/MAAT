function [OutLong,OutLat,TotRot]=convert_coo(Long,Lat,InCooType,OutCooType,JD,ObsCoo)
% Convert between different coordinates
% Package: celestial.coo
% Description: General coordinate convertor. Convert/precess
%              coordinate from/to Equatorial/galactic/Ecliptic
%              coordinates.
% Input  : - Matrix of longitudes [rad].
%          - Matrix of latitudes [rad].
%          - Type of input coordinates.
%            'J####.#' - equatorial, mean equinox of date
%            (default 'j2000.0'). [Julian years].
%            't####.#' - equatorial, true equinox of date
%            'tdate' - Equinox of provided JD.
%            'g','gal','galactic' - galactic.
%            'S' - Super galactic.
%            'c' - CMB dipole (see rotm_coo.m for details)
%            'e','ec','ecl','ecliptic' - Ecliptic with J2000.0 equinox.
%            'h','hor','azalt' - horizontal coordinates
%          - Type of outpt coordinates.
%            Like input. Default is 'g'.
%          - JD (for horizontal coordinates). Default is now.
%          - Observatory [Long,Lat] [rad] (for horizontal coordinates).
% Output : - Matrix of output longitude.
%          - Matrix of output latitude.
%          - Total rotation matrix.
% Tested : Matlab 2019b                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Lon,Lat,RM]=celestial.coo.convert_coo(rand(10,4),rand(10,4),'J2000.0','e')
%          [Lon,Lat]=celestial.coo.convert_coo(1.1,-0.2,'J2000.0','J2015.4')
% Reliable: 1
%--------------------------------------------------------------------------

RAD = 180./pi;

if (nargin<6)
    ObsCoo = [];
    if nargin<5
        JD = celestial.time.julday;
        if nargin<4
            OutCooType = 'g';
            if nargin<3
                InCooType  = 'J2000.0';
            end
        end
    end
end

switch lower(InCooType)
    case 'tdate'
        InCooType = sprintf('t%8.3f',convert.time(JD,'JD','J'));
end
switch lower(OutCooType)
    case 'tdate'
        OutCooType = sprintf('t%8.3f',convert.time(JD,'JD','J'));
end


LenInType  = length(InCooType);
LenOutType = length(OutCooType);

if ischar(InCooType)
   if (LenInType>1)
      InEquinox  = str2double(InCooType(2:LenInType));
      InEquinoxJD = 2451545.5 + 365.25.*(InEquinox - 2000);
      InCooType  = InCooType(1);
   end
else
    error('InCooType must be a string');
end

if ischar(OutCooType)
   if (LenOutType>1)
      OutEquinox = str2double(OutCooType(2:LenOutType));
      OutEquinoxJD = 2451545.5 + 365.25.*(OutEquinox - 2000);
      OutCooType  = OutCooType(1);
   end
else
    error('OutCooType must be a string');
end





RotM1 = diag([1 1 1]);

% calculate the first rotation matrix
switch lower(InCooType)
    case {'t','j'}
        if (InEquinox~=2000.0)
            % precess coordinates to J2000.0
            switch lower(InCooType)
                case {'j'}
                    % mean equinox ...
                    RotM1 = celestial.coo.rotm_coo('p',InEquinoxJD);
                case {'t'}
                    % true equinox ...
                    RotM1 = celestial.coo.rotm_coo('pd',InEquinoxJD);
                otherwise
                    error('Illegal InCooType');
            end
        end
    case {'g'}
        % convert to Equatorial J2000.0
        RotM1 = celestial.coo.rotm_coo('G',2451545.5);
    case {'s'}
        % convert to Equatorial J2000.0
        RotM1 = celestial.coo.rotm_coo('G',2451545.5)*rotm_coo('SGg');
    case {'c'}
        % convert to Equatorial
        error('not implemented')
    case {'e'}
        % convert to Equatorial J2000.0
        RotM1 = celestial.coo.rotm_coo('E',2451545.5);
    case {'h','a'}
        % convert to J2000
        if isempty(ObsCoo)
            error('when converting to/from AzAlt ObsCoo must be provided');
        end
        EqCoo = celestial.coo.horiz_coo([Long(:), Lat(:)],JD,ObsCoo,'e');
        Long  = EqCoo(:,1);
        Lat   = EqCoo(:,2);
    otherwise
        error('Unknown input coordinaytes type');
end

RotM2 = diag([1 1 1]);
% calculate the second rotation matrix
switch lower(OutCooType)
    case {'j','t'}
        if (OutEquinox~=2000.0)
            % precess coordinates from J2000.0
            switch lower(OutCooType)
                case {'j'}   
                    % mean equinox ...
                    RotM2 = celestial.coo.rotm_coo('P',OutEquinoxJD);
                case {'t'}   
                    % true equinox ...
                    RotM2 = celestial.coo.rotm_coo('Pd',OutEquinoxJD);
                otherwise
                    error('Illegal OutCooType');
            end           
        end
    case {'g'}
        % convert to galactic
        RotM2 = celestial.coo.rotm_coo('g',2451545.5);
     case {'s'}
        % convert to Super galactic
        RotM2 = celestial.coo.rotm_coo('gSG')*rotm_coo('g',2451545.5);
     case {'c'}
        % convert to CMB dipole
        error('not implemented')
     case {'e'}
        % convert to ecliptic
        RotM2 = celestial.coo.rotm_coo('e',2451545.5);
     case {'h','a'}
        % convert to horizontal
        if isempty(ObsCoo)
            error('when converting to/from AzAlt ObsCoo must be provided');
        end
        EqCoo = celestial.coo.horiz_coo([Long(:), Lat(:)],JD,ObsCoo(1:2),'h');
        Long  = EqCoo(:,1);
        Lat   = EqCoo(:,2);
     otherwise
        error('Unknown output coordinaytes type');
end

% convert coordinates to direction cosines
SizeLong = size(Long);
SizeLat  = size(Lat);
[CX,CY,CZ] = celestial.coo.coo2cosined(Long(:),Lat(:));
InCosDir = [CX(:), CY(:), CZ(:)];


% rotate coordinates
TotRot = RotM2*RotM1;
OutCosDir = TotRot*[InCosDir.'];

% convert coordinates from direction cosines
[OutLong,OutLat] = celestial.coo.cosined2coo(OutCosDir(1,:),OutCosDir(2,:),OutCosDir(3,:));

F0 = OutLong(:,1)<0;
OutLong(F0) = 2.*pi + OutLong(F0);

OutLong = reshape(OutLong,SizeLong);
OutLat  = reshape(OutLat,SizeLat);



