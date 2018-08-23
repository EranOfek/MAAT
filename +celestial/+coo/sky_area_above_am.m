function [SkyArea,LON]=sky_area_above_am(JD,Lat,AM,TimeVis)
% Calculate sky area observable during the night above a specific airmass.
% Package: celestial.coo
% Description: Calculate sky area observable during the night above
%              a specific airmass, and assuming each field is observable
%              for at least TimeVis hours.
% Input  : - JD.
%          - Latitute [rad].
%          - Airmass.
%          - Time visibility [hours].
% Output : - Sky area [deg^2]
%          - Length of night [fraction of day].
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Jul 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: for I=1:1:366, [SkyArea16(I),LON(I)]=celestial.coo.sky_area_above_am(VecJD(I),Lat,AM,TimeVis); end
% Reliable: 2
%--------------------------------------------------------------------------
import celestial.coo.*
import celestial.time.*
import celestial.SolarSys.*

RAD = 180./pi;

%Cat = tile_the_sky(60,30);
RA = 0;

Dec = (-30:1:90).'./RAD;

if (nargin==0)
   JD      = 2451545;
   Lat     = 33./RAD;
   AM      = 1.6;
   TimeVis = 1;  % [hour]
   
end


JD = floor(JD)+0.5;
[Time,~]=sun_rise_set(JD,[0 Lat],0,0);
UTend    = Time(2);
UTstart  = Time(8);
LON      = 1-(UTstart - UTend);  % length of night
LSTstart = lst(JD+UTstart,0).*2.*pi;  % [rad]
LSTend   = lst(JD+UTend,0).*2.*pi;    % [rad]

Alt = pi./2 - hardie_inv(AM);
% calculate southern declination limit
[Dec1,Dec2]=altha2dec(Alt,0.5.*TimeVis./24.*2.*pi,Lat);
DecLimitMin = min(Dec1,Dec2);
DecLimitMax = max(Dec1,Dec2);

HA  = alt2ha(Alt,Dec,Lat);
RAstart  = LSTstart - HA + TimeVis./24.*2.*pi;
RAend    = LSTend + HA - TimeVis./24.*2.*pi;

RAstart = (RAstart./(2.*pi) - floor(RAstart./(2.*pi))).*2.*pi;
RAend   = (RAend./(2.*pi) - floor(RAend./(2.*pi))).*2.*pi;

%Flag = RA>RAstart & RA<RAend & Dec>DecLimitMin & Dec<DecLimitMax;
RAdiff = RAend - RAstart;
    
FlagN = RAstart>RAend;
%Flag(FlagN) = (RA>RAstart(FlagN) | RA<RAend(FlagN)) & Dec>DecLimitMin & Dec<DecLimitMax;
RAdiff(FlagN) = (2.*pi-RAstart(FlagN)) + RAend(FlagN);

SkyArea = sum(RAdiff.*cos(Dec).*RAD);

