function [Long,Lat]=interp_diff_longlat(Time,LongLat,NewT,Deg,Check)
% Bessel interpolation of equally space time series of lon/lat coordinates
% Package: Util
% Description: Interpolate equally space time series of lon/lat coordinates
%              using 4th order Bessel differences interpolation.
% Input  : - 
%          - degree of differences, default is 4.
%          - Check if X is equally spaced {'y' | 'n'}, default is 'n'.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Long,Lat]=Util.interp.interp_diff_longlat(
% Reliable: 
%--------------------------------------------------------------------------


Def.Deg   = 4;
Def.Check = 'n';
if (nargin<6)
    Check = Def.Check;
    if (nargin<5)
        Deg = Def.Deg;
    end
end

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (size(LongLat,2)==2)
    [CX,CY,CZ] = celestial.coo.coo2cosined(LongLat(:,1),LongLat(:,2));
elseif (size(LongLat,2)==3)
    CX = LongLat(:,1);
    CY = LongLat(:,2);
    CZ = LongLat(:,3);
else
    error('unknown number of columns in LongLat - should be either 2 or 3');
end

[CXN] = Util.interp.interp_diff(Time,CX,NewT,Deg,Check);
[CYN] = Util.interp.interp_diff(Time,CY,NewT,Deg,Check);
[CZN] = Util.interp.interp_diff(Time,CZ,NewT,Deg,Check);

[Long,Lat] = celestial.coo.cosined2coo(CXN,CYN,CZN);



