function [CMg,CMi,Ds,De]=jup_meridian(JD)
% Low accuracy formula for Jupiter central meridian
% Package: celestial.SolarSys
% Description: Low accuracy formula for Jupiter central meridian.
% Input  : - Vector of JDs in TDT.
% Output : - matrix of central meridian (in deg.) [System_I, System_II]
%            The central meridian is calculated for the geometric disk.
%          - matrix of central meridian (in deg.) [System_I, System_II]
%            The central meridian is calculated for the illuminated disk.
%          - The planetocentric declination of the Sun. (deg.)
%          - The planetocentric declination of the Earth. (deg.)
% Reference: Meeus, J. 1991 in: Astronomical Algorithms.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Web example: http://astroclub.tau.ac.il/ephem/JovMap/index.php#JovMap
% Example: [CMg,CMi,Ds,De]=celestial.SolarSys.jup_meridian(2451545);
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;

d   = JD - 2451545.0;
V   = (172.74  + 0.00111588.*d)./RAD;
M   = (357.529 + 0.9856003.*d)./RAD;
N   = (20.020  + 0.0830853.*d + 0.329.*sin(V))./RAD;
J   = (66.115  + 0.9025179.*d - 0.329.*sin(V))./RAD;
A   = (1.915.*sin(M) + 0.020.*sin(2.*M))./RAD;
B   = (5.555.*sin(N) + 0.168.*sin(2.*N))./RAD;
K   = J + A - B;
R   = 1.00014 - 0.01671.*cos(M) - 0.00014.*cos(2.*M);
r   = 5.20872 - 0.25208.*cos(N) - 0.00611.*cos(2.*N);
Del = sqrt(r.^2 + R.^2 - 2.*r.*R.*cos(K));
Psi = asin(R.*sin(K)./Del);

% central meridian in deg.
SysI  = 210.98 + 877.8169088.*(d - Del./173) + (Psi - B).*RAD;
SysII = 187.23 + 870.1869088.*(d - Del./173) + (Psi - B).*RAD;

% reduce to [0-360]
SysI  = 360.*(SysI./360 - floor(SysI./360));
SysII = 360.*(SysII./360 - floor(SysII./360));

CMg   = [SysI, SysII];
CMi   = [CMg - sign(sin(K)).*57.3.*(sin(Psi./2)).^2];
%,...
%         CMg(:,2) -sign(sin(K)).*57.3.*(sin(Psi./2)).^2];

Lam   = (34.35 + 0.083091.*d + 0.329.*sin(V))./RAD + B;
Ds    = 3.12.*sin(Lam + 42.8./RAD);
De    = Ds - 2.22.*sin(Psi).*cos(Lam+22./RAD) - 1.30.*(r-Del).*sin(Lam-100.5./RAD)./Del;
