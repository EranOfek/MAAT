function [OutRA,OutDec]=coo_resolver(RA,varargin)
% Resolve coordinates or target name into RA/Dec 
% Package: celestial
% Description: Given coordinates (Lon/Lat) in any coordinate system or
%              format, or a target name convert the coordinates into
%              RA/Dec in some specific equinox and in deg/radians unuts.
% Input  : - Longitude in any format (e.g., deg, sexagesimal,...)
%            or system (e.g., equatorial, ecliptic,...).
%            If numeric than assume input is in degrees (unless 'InUnits'
%            is set to 'rad', for radians). If string and Latitute is not
%            empty, than assume input is in sexagesimal format
%            (e.g., 'HH:MM:SS.FFF' or 'HH MM SS.FFF'). If string and
%            Latitude is empty, or not given, than assume this is an object
%            name and try to use name resolver (e.g., NED/SIMBAD) to
%            convert the name to coordinates.
%            Available coordinate systems are listed in the 'InSys'
%            parameter. Default is J2000 equatorial.
%          - Latitude in any format or system.
%            If empty or not provided than assume Longitude is an opbject
%            name.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'InSys'    - Input coordinate system. Options are
%                         'j2000','el','gal' (see celestial.coo.coco for
%                         more options).
%                         Default is 'j2000'.
%            'InUnits'  - Input units in case of a numeric input
%                         coordinates {'deg'|'rad'}. Default is 'deg'.
%            'NameServer' - Name server function:
%                         @VO.name.server_ned [default].
%                         @VO.name.server_simbad
%            'OutUnits' - Output units {'deg'|'rad'}. Default is 'rad'.
% Output : - J2000.0 R.A.
%          - J2000.0 Dec.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: % RA/Dec [deg] to RA/Dec in radians
%          [OutRA,OutDec]=celestial.coo.coo_resolver(1,1);
%          % RA/Dec [rad] to RA/Dec in radians
%          [OutRA,OutDec]=celestial.coo.coo_resolver(1,1,'InUnits',rad');
%          % use NED name server (requires internet)
%          [OutRA,OutDec]=celestial.coo.coo_resolver('Deneb','NameServer',@VO.name.server_simbad)
%          [OutRA,OutDec]=celestial.coo.coo_resolver('M81',[])
%          [OutRA,OutDec]=celestial.coo.coo_resolver('M81','OutUnits','deg')
%          % input in sexagesimal
%          [OutRA,OutDec]=celestial.coo.coo_resolver('15:00:11.1','-13:01:11.1')
%          % input in galactic coordinates
%          [OutRA,OutDec]=celestial.coo.coo_resolver('15:00:11.1','-13:01:11.1','InSys','gal')
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin.*0.5)==floor(nargin.*0.5)
    % even number of arguments
    % RA, Dec are provided
    Dec = varargin{1};
    varargin = varargin(2:end);
else
    % odd number of arguments
    Dec = [];
end




DefV.InSys                = 'j2000';     % 'j2000' | 'gal' | 'ec'
DefV.InUnits              = 'deg';    % 'rad'
DefV.NameServer           = @VO.name.server_ned;    % @VO.name.server_simbad 
DefV.OutUnits             = 'rad';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% treat empty Dec - i.e., object name to RA/Dec
if isempty(Dec)
    if (ischar(RA))
        Name = RA;
        [RA, Dec] = InPar.NameServer(Name);
        % output is in deg [J2000]
        InPar.InUnits   = 'deg';
        InPar.InSys     = 'j2000.0';
        InPar.InEquinox = 'j2000';
        
    else
        error('Dec is not provided so the first input argument is expected to be a string');
    end
end

if (isnumeric(RA) && isnumeric(Dec))
    % convert input coordinates to radians
    Factor = convert.angular(InPar.InUnits,'rad');
    RA     = RA.*Factor;
    Dec    = Dec.*Factor;
end

% treat Coordinates format
RA     = celestial.coo.convertdms(RA,'gH','r');
Dec    = celestial.coo.convertdms(Dec,'gD','R');

% treat coordinate system
Coo = celestial.coo.coco([RA, Dec], InPar.InSys, 'j2000');

% convert to outoput units
Factor = convert.angular('rad',InPar.OutUnits);
Coo    = Coo.*Factor;
OutRA  = Coo(:,1);
OutDec = Coo(:,2);


    
