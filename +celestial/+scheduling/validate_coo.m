function Flag=validate_coo(HA,Dec,varargin)
% Validate HA, Dec within Az,Alt and HA,Dec ranges.
% Package: +celestial.scheduling
% Description: Given a vectors of HA/Dec check for each coordinate if it
%              complies with Az,Alt limit and withing HA/Dec ranges.
% Input  : - Matrix of HA [rad].
%          - Matrix of Dec [rad].
%          * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'Lon' - Observatory East longitude [rad]. Default is 35./RAD.
%            'Lat' - Observatory North Latitude [rad]. Default is 30./RAD.
%            'AzAltLimit' - A two column matrix of [Az Alt] limits [deg].
%                   First line must be for Az=0 and last line for Az=360.
%                   Default is [0 15; 90 15; 180 15; 270 15; 360 15]
%            'RangeHA' - Range of valid HA (-180 to 180) [deg].
%                   Default is [-120 120].
%            'RangeDec' - Range of valid Dec [deg].
%                   Default is [-30 90].
%            'InterpMethod' - Table interpolation method.
%                   Default is 'linear'.
% Output : - A matrix of logicals indicating if target is in range (true).
% Example: Flag=celestial.scheduling.validate_coo(HA,Dec)


RAD = 180./pi;

InPar = inputParser;
addOptional(InPar,'Lon',35./RAD);
addOptional(InPar,'Lat',30./RAD);
addOptional(InPar,'AzAltLimit',[0 15; 90 15; 180 15; 270 15; 360 15]);
addOptional(InPar,'RangeHA',[-120 120]);  % [deg] in range -pi to pi
addOptional(InPar,'RangeDec',[-30 90]); % [deg]
addOptional(InPar,'InterpMethod','linear');
parse(InPar,varargin{:});
InPar = InPar.Results;

if InPar.AzAltLimit(1,1)~=0 && InPar.AzAltLimit(end,1)~=360
    error('AzAltLimit must start with Az 0 and end at Az 360');
end

% convert HA to range -pi to pi
HA = mod(HA,2.*pi);
HA(HA>pi) = HA(HA>pi) - 2.*pi;

[Az,Alt] = celestial.coo.hadec2azalt(HA,Dec,InPar.Lat);

AltLimit = interp1(InPar.AzAltLimit(:,1),InPar.AzAltLimit(:,2),Az,InPar.InterpMethod);

Flag = Alt>(AltLimit./RAD) & HA>=(InPar.RangeHA(1)./RAD) & HA<=(InPar.RangeHA(2)./RAD) & Dec>=(InPar.RangeDec(1)./RAD) & Dec<=(InPar.RangeDec(2)./RAD);




