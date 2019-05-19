function Sun=get_sun(JD,GeodCoo)
% Get Sun position (low accuracy)
% Package: celestial.SolarSys
% Description: Get Sun position (low accuracy).
% Input  : - JD, or date (see jd.m for available formats).
%          - Geodetic coordinates [Long, Lat] in radians.
% Output : - Structure containing Sun position, with the following fields:
%            .RA     - RA [radians]
%            .Dec    - Dec [radians]
%            .Az     - Azimuth [radians]
%            .Alt    - Altitude [radians]
%            .dAzdt  - dAz/dt [radians/s] 
%            .dAltdt - dAlt/dt [radians/s] 
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Jun 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sun=celestial.SolarSys.get_sun(2451545,[1 1]);
% Reliable: 1
%--------------------------------------------------------------------------
DT = 1./86400;   % 1 s [days]
if (size(JD,2)==1)
   % already in JD
else
  JD = celestial.time.julday(JD).';
end

[Sun.RA, Sun.Dec]   = celestial.SolarSys.suncoo(JD,'a');
[SunN.RA, SunN.Dec] = celestial.SolarSys.suncoo(JD+DT,'a');
I = find(Sun.RA<0);
if (isempty(I)==1),
   % do nothing
else
   Sun.RA(I) = 2.*pi + Sun.RA(I);
end
HorizCoo  = celestial.coo.horiz_coo([Sun.RA, Sun.Dec],JD,GeodCoo,'h');
HorizCooN = celestial.coo.horiz_coo([SunN.RA, SunN.Dec],JD+DT,GeodCoo,'h');
Sun.Az    = HorizCoo(:,1);
Sun.Alt   = HorizCoo(:,2);

% derivatives in [Rad/sec]
Sun.dAzdt  = (HorizCooN(:,1) - HorizCoo(:,1));    
Sun.dAltdt = (HorizCooN(:,2) - HorizCoo(:,2));
