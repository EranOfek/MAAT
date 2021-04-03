function [Flag,FlagRes,Data]=is_coordinate_ok(RA,Dec,JD,varargin)
% Check that coordinates satisfy some observability conditions
% Package: celestial.coo
% Description: Check that J2000 equatorial coordinates satisfy some
%              observability conditions including Az, Alt, HA.
% Input  : - J2000.0 R.A. [radians].
%          - J2000.0 Dec. [radians].
%          - JD. If emepty than use current JD.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Lon'
%            '
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Flag,FlagRes]=celestial.coo.is_coordinate_ok(1,0)
% Reliable: 
%--------------------------------------------------------------------------


if (nargin==2)
    JD = [];
end

RAD = 180./pi;

DefV.Lon                  = 35./RAD;
DefV.Lat                  = 32./RAD;
DefV.AzAltConst           = [0 15; 90 15; 270 15; 360 15]./RAD;
DefV.HAConst              = [-0.75.*pi,  0.75.*pi];
DefV.AltMinConst          = 15./RAD;
DefV.AltMaxConst          = 90./RAD;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (isempty(JD))
    JD = celestial.time.julday;
end

LST = celestial.time.lst(JD,InPar.Lon,'a');  % [frac of day]
LST = LST.*2.*pi;  % [rad]

% calc HA
HA = LST - RA;  % [rad]
% calc Az/Alt
% do we need to convert RA/Dec to equinox of date?
AzAlt = celestial.coo.horiz_coo([RA, Dec],JD,[InPar.Lon, InPar.Lat],'h');
Az    = AzAlt(:,1);
Alt   = AzAlt(:,2);
% calc AM
AM = celestial.coo.hardie(pi./2-Alt);

% if (Az<0 || Az>(2.*pi))
%     error('Illegal Az')
% end

% set HA range to -pi to pi:
HA = mod(HA,2.*pi);
FF = HA > pi;
HA(FF) = HA(FF) - 2.*pi;


% if (HA>pi || HA<-pi)
%     error('Illgal HA');
% end

Data.LST = LST;
Data.HA  = HA;
Data.Az  = Az;
Data.Alt = Alt;
Data.AM  = AM;

FlagRes.AzAlt = Alt>interp1(InPar.AzAltConst(:,1),InPar.AzAltConst(:,2),Az);
FlagRes.HA    = HA>InPar.HAConst(1) & HA<InPar.HAConst(2);
FlagRes.Alt   = Alt>InPar.AltMinConst & Alt<=InPar.AltMaxConst;

Flag = FlagRes.AzAlt & ...
       FlagRes.HA    & ...
       FlagRes.Alt;
