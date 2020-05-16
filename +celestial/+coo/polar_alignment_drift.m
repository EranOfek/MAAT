function [dRt_dt,dDt_dt]=polar_alignment_drift(H,D,Ht,Dt,Phi,Psi,Beta)
% SHORT DESCRIPTION HERE
% Package: celestial
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [dHt_dt,dDt_dt]=celestial.coo.polar_alignment_drift(1,0,1,0,32./RAD,0./RAD,1./RAD)
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;


% D - true Dec
% Dt - telescope Dec.
% R - true RA
% Rt - telescope RA
% H - true HA
% Ht - telescope HA
% Phi - true altitude of NCP (latitude)
% Phit - telescope polar altitude
% Z - zenith distance of star
% Beta - distance between NCP and telescope pole
% Psi - HA of telescope pole

Phit = asin(sin(Phi)*cos(Beta) + cos(Phi).*sin(Beta).*cos(Psi));
Alt = celestial.coo.ha2alt(H,D,Phi);
Z   = pi./2 - Alt;

dH_dt = 15; % "/s
dH_dt = dH_dt./(3600.*RAD);

dDt_dt = dH_dt .* (cos(D).*sin(Beta).*sin(Psi-H))./cos(Dt);

Ht = acos(cos(Z)./(cos(Dt).*cos(Phit)) - tan(Dt).*tan(Phit));

% note that there is an error in the last term of Equation 4 in Markworth
dHt_dt = dH_dt .* cos(D).*cos(Phi).*sin(H)./(cos(Dt).*cos(Phit).*sin(Ht)) - sec(Dt).^2.*dDt_dt.*cos(Z).*sin(Dt)./(sin(Ht).*cos(Phit)) + ...
         dDt_dt.*tan(Phit)./(cos(Dt).^2 .* sin(Ht));

     
dDt_dt = dDt_dt.*3600.*RAD;
dRt_dt = (dHt_dt - dH_dt).*3600.*RAD;
