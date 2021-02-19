function [DistRA,DistDec,Aux]=convert2equatorial(Long,Lat,varargin)
% Convert coordinates/name to apparent equatorial coordinates. 
% Package: celestial
% Description: Given a coordinates in some coordinate system or equinox,
%              or an object name, convert it to euatorial coordinates that
%              includes the atmospheric refraction correction and optional
%              telescope distortion model (T-point model).
% Input  : - Longitude in some coordinate system, or object name.
%            Longitude can be either sexagesimal coordinates or numeric
%            calue in degress (or radians if InputUnits='rad').
%            Object name is converted to coordinates using either SIMBAD,
%            NED or JPL horizons.
%          - Like the first input argument, but for the latitude.
%            If empty, or not provided, than the first argument is assumed
%            to be an object name.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            However, if only one parameter is provided, it is treated as
%            the value of 'InCooType'.
%            'InCooType'  - Input coordinates frame:
%                           'a' - Az. Alt.
%                           'g' - Galactic.
%                           'e' - Ecliptic
%                           'ha' - Hour Angle/Declination at true equinox
%                           of date.
%                           'jha' - Hour Angle/Declination at J2000.0.
%                           - A string start with J (e.g., 'J2000.0').
%                           Equatorial coordinates with mean equinox of
%                           date, where the year is in Julian years.
%                           -  A string start with t (e.g., 't2020.5').
%                           Equatorial coordinates with true equinox of
%                           date.
%                           Default is 'J2000.0'
%            'OutCooType' - 'tdate' | 'J2000.0'. Default is 'J2000.0'.
%            'NameServer' - ['simbad'] | 'ned' | 'jpl'.
%            'JD'         - Julian day. This is used for horizontal
%                           coordinates and H.A.
%                           Default is now (i.e., celestial.time.julday)
%            'ObsCoo'     - Observer Geodetic position.
%                           [East Long (deg), Lat (deg), Height (meters)]
%            'HorizonsObsCode' - In case solar system ephemerides is
%                           requested, then this is the Horizons observatory
%                           code. Default is '500' (geocentric observer).
%            'DistFun'    - Distortion function handle.
%                           The function is of the form:
%                           [DistHA,DistDec]=@Fun(HA,Dec), where all the
%                           input and output are in degrees.
%                           Default is empty. If not given return [0,0].
%            'InputUnits' - Default is 'deg'.
%            'OutputUnits'- Default is 'deg'
%            'ApplyRefraction' - Default is true.
%            'Temp'       - Default is 15 C.
%            'Wave'       - Default is 5500 Ang.
%            'PressureHg' - Default is 760 mm Hg.
% Output : - Apparent R.A. (or HA, if InCooType is 'ha', or 'jha')
%          - Apparent Dec.
%          - A structure containing the intermidiate values.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [DistRA,DistDec,Aux]=celestial.coo.convert2equatorial(1,1)
%          celestial.coo.convert2equatorial('12:00:00','+20:00:00');
%          celestial.coo.convert2equatorial('M31');
%          celestial.coo.convert2equatorial('9804;',[],'NameServer','jpl')
%          celestial.coo.convert2equatorial('12:00:00','+20:00:00','InCooType','J2000.0');
%          celestial.coo.convert2equatorial('12:00:00','+20:00:00','J2000.0');
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;

if nargin.*0.5~=floor(nargin.*0.5)
    varargin = {'InCooType',varargin{:}};
end

%OutCooType = 'J2000.0';  % or 'tdate'
if nargin<2
    Lat = [];
end

InPar = inputParser;


addOptional(InPar,'InCooType','J2000.0');   % 'eq' | 'gal' | 'ecl' | 'horizon'
addOptional(InPar,'OutCooType','J2000.0');   % 'eq' | 'gal' | 'ecl' | 'horizon'
addOptional(InPar,'NameServer','simbad');  % 'simbad' | 'ned' | 'jpl'
addOptional(InPar,'JD',celestial.time.julday);  % time for solar system ephemerids
addOptional(InPar,'ObsCoo',[35 30.6 800]);  % 
addOptional(InPar,'HorizonsObsCode','500');  % 500 geocentric

addOptional(InPar,'InputUnits','deg');  
addOptional(InPar,'OutputUnits','deg');  
addOptional(InPar,'DistFun',[]);  
addOptional(InPar,'ApplyRefraction',true);
addOptional(InPar,'Temp',15);  % C
addOptional(InPar,'Wave',5500);  % Ang
addOptional(InPar,'PressureHg',760);  % mm Hg

parse(InPar,varargin{:});

InPar = InPar.Results;

InPar.ObsCoo(1:2) = convert.angular('deg','rad',InPar.ObsCoo(1:2));




if isempty(Lat)
    % assume user supplied object name
    switch lower(InPar.NameServer)
        case 'simbad'
            [Long,Lat] = VO.name.server_simbad(Long);  % output in degrees
            InPar.InputUnits = 'deg';
        case 'ned'
            [Long,Lat] = VO.name.server_ned(Long);  % output in degrees
            InPar.InputUnits = 'deg';
        case 'jpl'
            % % semicolumn telss horizon its a small body - e.g., '499;'
            [JCat]=celestial.SolarSys.jpl_horizons('ObjectInd',Long,'StartJD',InPar.JD-2,'StopJD',InPar.JD+2,'StepSizeUnits','h','CENTER',InPar.HorizonsObsCode);
            [Long,Lat] = Util.interp.interp_diff_longlat(JCat.Cat(:,JCat.Col.JD),...
                                                [JCat.Cat(:,JCat.Col.RA),JCat.Cat(:,JCat.Col.Dec)],...
                                                InPar.JD);     
            InPar.InputUnits = 'rad';
        otherwise
            error('Unknown NameServer option');
    end
    
    % assume input is J2000.0
    InPar.InCooType     = 'J2000.0';
    
    % convert coordinates to radians
    Long = convert.angular(InPar.InputUnits,'rad',Long);
    Lat  = convert.angular(InPar.InputUnits,'rad',Lat);
    
else
    if ischar(Long)
        Long = celestial.coo.convertdms(Long,'SH','r');
    else
        Long = convert.angular(InPar.InputUnits,'rad',Long);
    end
    if ischar(Lat)
        Lat = celestial.coo.convertdms(Lat,'SD','r');
    else
        Lat  = convert.angular(InPar.InputUnits,'rad',Lat);
    end
    
end

% calculate LST
LST = celestial.time.lst(InPar.JD,InPar.ObsCoo(1),'a');  % fraction of day

switch lower(InPar.InCooType)
    case 'ha'
        % HA/Dec at true equinox of date
        JulianYear = convert.time(InPar.JD,'JD','J');
        InPar.InCooType = sprintf('J%7.3f',JulianYear);
        
        % HA = LST - RA
        % in this case Long input is HA, and we convert it to RA
        Long = 2.*pi.*LST - Long;  % Long is now RA
        
    case 'jha'
        % HA/Dec at J2000.0
        InPar.InCooType = 'J2000.0';
        
        % HA = LST - RA
        % in this case Long input is HA, and we convert it to RA
        Long = 2.*pi.*LST - Long;  % Long is now RA
        
    otherwise
        % do nothing
end


% convert input coordinates to RA, Dec in output equinox
[TrueRA,TrueDec] = celestial.coo.convert_coo(Long,Lat,InPar.InCooType,InPar.OutCooType,InPar.JD,InPar.ObsCoo);
TrueHA = 2.*pi.*LST - TrueRA;

% applay atmospheric refraction
% in order for this step to be exact: the HA/Dec should be in equinox of
% date
[TrueAz,TrueAlt]=celestial.coo.hadec2azalt(TrueHA,TrueDec,InPar.ObsCoo(2));

%[TrueAz,TrueAlt] = celestial.coo.convert_coo(TrueRA,TrueDec,InPar.OutCooType,'azalt',InPar.JD,InPar.ObsCoo);
[Refraction]     = celestial.coo.refraction_wave(TrueAlt,InPar.Wave,InPar.Temp,InPar.PressureHg);
% add refraction to TrueAlt
AppAz  = TrueAz;
if InPar.ApplyRefraction
    AppAlt = TrueAlt + Refraction;
else
    AppAlt = TrueAlt;
end
% return to equatorial coordinates
%[AppRA,AppDec] = celestial.coo.convert_coo(AppAz,AppAlt,'azalt',InPar.OutCooType,InPar.JD,InPar.ObsCoo);
[AppHA,AppDec] = celestial.coo.azalt2hadec(AppAz,AppAlt,InPar.ObsCoo(2));
AppRA = 2.*pi.*LST - AppHA;  % [rad]

% applay distortions
% calculate HA = LST - RA
%AppHA = 2.*pi.*LST - AppRA;   % [rad]
% call distortions function
if isempty(InPar.DistFun)
    DeltaDistHA  = 0; % deg
    DeltaDistDec = 0; % deg
    
    DistHA  = AppHA  + DeltaDistHA./RAD;
    DistRA  = AppRA  + DeltaDistHA./RAD;
    DistDec = AppDec + DeltaDistDec./RAD;

else
    [DistHA, DistDec] = InPar.DistFun(AppHA,AppDec);
end

Aux.JD        = InPar.JD;
Aux.LST       = LST;  % fraction of day
Aux.TrueRA    = TrueRA;  % rad
Aux.TrueDec   = TrueDec; % rad
Aux.AppRA     = AppRA;
Aux.AppDec    = AppDec;
Aux.AppHA     = AppHA;
Aux.DistHA    = DistHA;  % rad
Aux.DistRA    = DistRA;  % rad
Aux.DistDec   = DistDec; % rad
Aux.TrueAz    = TrueAz;
Aux.TrueAlt   = TrueAlt;
Aux.AppAz     = AppAz;
Aux.AppAlt    = AppAlt;
Aux.DeltaDistHA  = DeltaDistHA./RAD;
Aux.DeltaDistDec = DeltaDistDec./RAD;
Aux.Refraction   = Refraction;

% convert RA back to HA if needed
switch lower(InPar.InCooType)
    case {'ha','jha'}
        % HA = LST - RA
        DistRA = 2.*pi.*LST - DistRA;  % DistRA is now HA
    otherwise
        % do nothing
end        
        

% convert back to output units
DistRA  = convert.angular('rad',InPar.OutputUnits,DistRA);
DistDec = convert.angular('rad',InPar.OutputUnits,DistDec);







