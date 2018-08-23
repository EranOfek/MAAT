function Out=get_horizons(StartJD,EndJD,Object,varargin)
% Get an ephemerides for a solar system body from JPL horizons
%--------------------------------------------------------------------------
% get_horizons function                                              ephem
% Description: get an ephemerides for a solar system body from the JPL
%              horizons system.
%              OBSOLETE: Use jpl_horizons.m instead
% Input  : - Start JD or date [D M Y Frac], [D M Y H M S], [D M Y].
%            Default time scale is UTC. Use 'TimeScale' argument for
%            other time scales (e.g., 'TT').
%          - End JD or date.
%          - Object name:
%            'Mercury', 'Venus', 'Mars',
%            'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'
%            'Moon','Sun'
%            Alternatively, a vector of orbital elements
%            [e, q, T, Long_AscNode, Arg_Peri, Inc] (au, JD, deg)
%            or
%            ...
%            
%          * Arbitrary number of pairs of input arguments:
%            ...,'keyword','value',...
%            The following keywords are available:
%            'Geod'      - Geodetic coordinates.
%                          If 'Geocentric' (default), the calculate
%                          geocentric coordinates.
%                          Otherwise [Long, Lat, Height].
%                          Where Long and Lat are in deg, and Height in km.
%                          Coordinates are given relative to WGS84.
%            'StepSize'  - Tabulation step size in days, default is 1 day.
%            'StepUnit'  - {'d' | 'h' | 'm'}, default is 'd'.
%            'TimeScale' - Time scale {'UT' | 'TT'}, default is 'UT'.
%            'Type'      - planet name type {'mp' | 'p'}, default is 'mp'.
%
%
% Output : - Structure containing the object ephemerides:
%            .JD        - Julian day
%            .J_RA      - J2000.0 astrometric R.A. [rad]
%            .J_Dec     - J2000.0 astrometric Dec. [rad]
%            .J_RA_sex  - J2000.0 astrometric R.A. [sex]
%            .J_Dec_sex - J2000.0 astrometric Dec [sex]
%            .a_RA      - Apparent R.A. with respect to true Equinox of date [rad].
%                         Corrected for light-time, the gravitational deflection
%                         of light, stellar aberration, precession and nutation.
%            .a_Dec     - Apparent Dec. with respect to true Equinox of date [rad].
%            .a_RA_sex  - Apparent R.A. with respect to true Equinox of date [sex].
%            .a_Dec_sex - Apparent Dec. with respect to true Equinox of date [sex].
%            .Mag       - Apparent V-band magnitude.
%            .SB        - Surface magnitude [mag/arcsec^2]
%            .Illum     - Illuminated fraction
%            .AngDiam   - Angular diamater ["]
%            .NP_PA     - PA of north pole [rad].
%            .NP_Dist   - Ang. distance to north pole ["].
%                         Negative value indicate NP on hidden hemisphere.
%            .r         - Target-Sun distance [au]
%            .Delta     - Target-Observer distance [au].
%            .AngSOT    - Sun-Observer-Target angle [rad]
%                         Apparent solar elongation, if negative the target center
%                         behind the Sun.
%            .AngSOTd   - +1 tragets leads the Sun (morning sky)
%                         -1 tragets trail the Sun (evening sky)
%            .AngSTO    - Sun-Target-Observer abgle [rad]
%                         Apparent phase angle of target.
%            .AngTOM    - Target-Observer-Moon angle [rad]
%            .MoonIllum - Moon illuminated fraction
%            .Const     - Constellation
% Reference: http://ssd.jpl.nasa.gov/horizons_batch.cgi
%            http://ssd.jpl.nasa.gov/?horizons_doc#specific_quantities
%            ftp://ssd.jpl.nasa.gov/pub/ssd/horizons_batch_example.long
% Tested : Matlab 7.3
%     By : Eran O. Ofek                   Jul 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Data=celestial.SolarSys.get_horizons([1 1 2000],[1 2 2000],'Mars');
%          Data=celestial.SolarSys.get_horizons([1 1 2000],[1 2 2000],'500%3B','Type','p');
% Reliable: 2
%--------------------------------------------------------------------------
RAD = 180./pi;

if (length(StartJD)==1)
   % given in JD
   %StartDate = jd2date(floor(StartJD)+0.5);
   StartDate = celestial.time.jd2date(StartJD,'h');
else
   if (length(StartJD)==3)
      StartDate = [StartJD, 0 ,0];
   elseif (length(StartJD)==5)
      StartDate = [StartJD];
   elseif (length(StartJD)==6)
      StartDate = StartDate(1:5);
   else
      error('Illegal date format');
   end
end
if (length(EndJD)==1)
   % given in JD
   %EndDate = jd2date(floor(EndJD)+0.5);
   EndDate = jd2date(EndJD,'h');
else
   if (length(EndJD)==3)
      EndDate = [EndJD, 0 ,0];
   elseif (length(EndJD)==5)
      EndDate = EndJD;
   elseif (length(EndJD)==6)
      EndDate = EndDate(1:5);
   else
      error('Illegal date format');
   end
end

if (StartDate(3)<1000 && StartDate(3)>0),
   StartDateStr = sprintf('%04dAD-%d-%d%%%20%02d:%02d',StartDate([3 2 1 4 5]));
elseif (StartDate(3)<=0),
   StartDateStr = sprintf('%04dBC-%d-%d%%20%02d:%02d',StartDate([3 2 1 4 5]));
else
   StartDateStr = sprintf('%04d-%d-%d%%20%02d:%02d',StartDate([3 2 1 4 5]));
end

if (EndDate(3)<1000 && EndDate(3)>0),
   EndDateStr   = sprintf('%04dAD-%d-%d%%20%02d:%02d',EndDate([3 2 1 4 5]));
elseif (EndDate(3)<=0),
   EndDateStr   = sprintf('%04dBC-%d-%d%%20%02d:%02d',EndDate([3 2 1 4 5]));
else
   EndDateStr   = sprintf('%04d-%d-%d%%20%02d%:02d',EndDate([3 2 1 4 5]));
end

%StartDateStr='5BC-01-01'
%EndDateStr='1BC-01-01'

%--- set default values ---
Geod      = 'Geocentric';
StepSize  = 1;
StepUnit  = 'd';
TimeScale = 'UT';
Type      = 'mp';
Narg = length(varargin);
if (Narg./2~=floor(Narg./2)),
   error('Illegal number of input arguments');
end
for Iarg=1:2:Narg-1,
   switch lower(varargin{Iarg})
    case 'stepsize'
       StepSize   = varargin{Iarg+1};
    case 'stepunit'
       StepUnit   = varargin{Iarg+1};
    case 'timescale'
       TimeScale  = varargin{Iarg+1};
    case 'geod'
       Geod       = varargin{Iarg+1};
    case 'type'
       Type       = varargin{Iarg+1};
    otherwise
       error('Unknown input argument option');
   end
end

if (ischar(Geod)==1),
   switch lower(Geod)
    case 'geocentric'
       Topo = 'n';
    otherwise
        error('Unknown Geod option');
   end
else
   Topo = 'y';
end

PlType = 'p';
switch lower(Object)
 case 'mercury'
    PlanetCode = 199;
 case 'venus'
    PlanetCode = 299;
 case 'mars'
    PlanetCode = 499;
 case 'jupiter'
    PlanetCode = 599;
 case 'saturn'
    PlanetCode = 699;
 case 'uranus'
    PlanetCode = 799;
 case 'neptune'
    PlanetCode = 899;
 case 'pluto'
    PlanetCode = 999;
 case 'moon'
    PlanetCode = 301;
 case 'sun'
    PlanetCode = 10;
 case 'io'
    PlanetCode = 501;
 case 'europa'
    PlanetCode = 502;
 case 'ganymede'
    PlanetCode = 503;
 case 'callisto'
    PlanetCode = 504;
 case 'titan'
    PlanetCode = 606;
 case {'ceres'}
    PlanetCode = '1%3B';
 case {'juno'}
    PlanetCode = '2%3B';
 case {'pallas'}
    PlanetCode = '3%3B';
 case {'vesta'}
    PlanetCode = '4%3B';

 otherwise
    % other object - use as planet code
    PlanetCode = Object;
    PlType = 'p';
    PlType = 'mp';
    PlType = Type;
end

if (ischar(PlanetCode)==1),
   PlanetCode = strrep(PlanetCode,' ','%20');
end

Quant      = '1,2,9,10,13,17,19,20,23,24,25,29';

switch Topo
 case 'n'
    if (ischar(PlanetCode)==0)
       URL=sprintf('http://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND="%d"&MAKE_EPHEM="YES"&TABLE_TYPE="OBSERVER"&START_TIME="%s%s"&STOP_TIME="%s"&STEP_SIZE="%d%%20%s"&QUANTITIES="%s"&CSV_FORMAT="YES"',PlanetCode,StartDateStr,TimeScale,EndDateStr,StepSize,StepUnit,Quant);
    else
       URL=sprintf('http://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND="%s"&MAKE_EPHEM="YES"&TABLE_TYPE="OBSERVER"&START_TIME="%s%s"&STOP_TIME="%s"&STEP_SIZE="%d%%20%s"&QUANTITIES="%s"&CSV_FORMAT="YES"',PlanetCode,StartDateStr,TimeScale,EndDateStr,StepSize,StepUnit,Quant);
    end
 case 'y'
    if (ischar(PlanetCode)==0)
       URL=sprintf('http://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND="%d"&MAKE_EPHEM="YES"&TABLE_TYPE="OBSERVER"&START_TIME="%s%s"&STOP_TIME="%s"&STEP_SIZE="%d%%20%s"&QUANTITIES="%s"&CSV_FORMAT="YES"&CENTER="COORD"&COORD_TYPE="GEODETIC"&SITE_COORD="%09.5f,%09.5f,%09.4f"',PlanetCode,StartDateStr,TimeScale,EndDateStr,StepSize,StepUnit,Quant,Geod(1), Geod(2), Geod(3));
    else
       URL=sprintf('http://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND="%s"&MAKE_EPHEM="YES"&TABLE_TYPE="OBSERVER"&START_TIME="%s%s"&STOP_TIME="%s"&STEP_SIZE="%d%%20%s"&QUANTITIES="%s"&CSV_FORMAT="YES"&CENTER="COORD"&COORD_TYPE="GEODETIC"&SITE_COORD="%09.5f,%09.5f,%09.4f"',PlanetCode,StartDateStr,TimeScale,EndDateStr,StepSize,StepUnit,Quant,Geod(1), Geod(2), Geod(3));
    end
 otherwise
    error('Unknown Topo option');
end

%URL

Str=urlread(URL);
Str=strrep(Str,'n.a.','NaN');


Is = strfind(Str,'$$SOE');
Ie = strfind(Str,'$$EOE');

ClearStr = Str(Is+7:Ie-1);
FID = fopen('try.horizon','w');
fprintf(FID,'%s',ClearStr);
fclose(FID);


%                 Date emp  emp  RA   Dec  RA   Dec  Mag, SB Ill  AD,  PA,  NPd
% r rdot delta deltadot, S-O-T, trailing/leading, S-T-O, T-O-M, Moon_Illum, Const


switch lower(PlType)
 case 'p' 
    Data = textscan(ClearStr,'%s %s %s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %s %f %f %f %s\n','delimiter',',');
 case 'mp'
    %                         1  2  3  4  5  6   7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
    Data = textscan(ClearStr,'%s %s %s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %s %f %f %f %s\n','delimiter',',');
%    Data = textscan(ClearStr,'%s %s %s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %s %f %f %f %s %s\n','delimiter',',');
 otherwise
    error('Unknown Object Type');
end

[Date,JD]=Util.string.read_date(Data{1},...
                    {1,4,'%d'; 6,8,'%s'; 10, 11, '%d'; 13, 14, '%d'; 16,17,'%d'},...
                    {'Year','MonthNameShort','Day','Hour','Min'});

Data{1} = JD;

Out.JD        = JD.';
Out.J_RA      = celestial.coo.convertdms(Data{4},'SHb','r');
Out.J_Dec     = celestial.coo.convertdms(Data{5},'SDb','R');
Out.J_RA_sex  = Data{4};
Out.J_Dec_sex = Data{5};
Out.a_RA      = celestial.coo.convertdms(Data{6},'SHb','r');
Out.a_Dec     = celestial.coo.convertdms(Data{7},'SDb','R');
Out.a_RA_sex  = Data{6};
Out.a_Dec_sex = Data{7};
Out.Mag       = Data{8};

switch lower(PlType)
 case 'p'
    Out.SB        = Data{9};
    Out.Illum     = Data{10}./100;
    Out.AngDiam   = Data{11};        % ["]
    Out.NP_PA     = Data{12}./RAD;   % [rad]   % relative to the true celestial north pole of date
    Out.NP_Dist   = Data{13};        % ["]     % negative value indicate NP on hidden hemisphere
    Out.r         = Data{14};        % [au]
    Out.Delta     = Data{16};        % [au]
    Out.AngSOT    = Data{18}./RAD;   % [rad]   % Sun-Observer-Target angle
                                     % apparent solar elongation
                                     % if negative the target center behind the Sun
    Out.AngSOTd   = Util.cell.isempty_cell(strfind(Data{19},'/T')).*2-1;
                                     % +1 tragets leads the Sun (morning sky)
                                     % -1 tragets trail the Sun (evening sky)
    Out.AngSTO    = Data{20}./RAD;   % [rad]   % Sun-Target-Observer abgle
                                     % apparent phase angle of target
    Out.AngTOM    = Data{21}./RAD;   % [rad]  Target-Observer-Moon angle
    Out.MoonIllum = Data{22}./100;   % Moon illuminated fraction
    Out.Const     = Data{23};        % Constellation

 case 'mp'
    Out.SB        = NaN;
    Out.Illum     = Data{9}./100;
    Out.AngDiam   = Data{10};        % ["]
    Out.NP_PA     = Data{11}./RAD;   % [rad]   % relative to the true celestial north pole of date
    Out.NP_Dist   = NaN;
    Out.r         = Data{13};        % [au]
    Out.Delta     = Data{15};        % [au]
    Out.AngSOT    = Data{17}./RAD;   % [rad]   % Sun-Observer-Target angle
                                     % apparent solar elongation
                                     % if negative the target center behind the Sun
    Out.AngSOTd   = Util.cell.isempty_cell(strfind(Data{18},'/T')).*2-1;
                                     % +1 tragets leads the Sun (morning sky)
                                     % -1 tragets trail the Sun (evening sky)
    Out.AngSTO    = Data{19}./RAD;   % [rad]   % Sun-Target-Observer abgle
                                     % apparent phase angle of target

    Out.AngTOM    = Data{20}./RAD;   % [rad]  Target-Observer-Moon angle
    Out.MoonIllum = Data{21}./100;   % Moon illuminated fraction
    Out.Const     = Data{22};        % Constellation

 otherwise
    error('Unknown Object Type');
end
