function [dHt_dt,dDt_dt,dRt_dt,Ht,Dt]=polar_alignment_drift(H,D,Phi,Psi,Beta)
% Calculate the RA/Dec drift due to equatorial polar alignemnt error.
% Package: celestial
% Description: 
% Input  : - HA of target [rad].
%          - Dec of target [rad].
%          - True altitude of NCP (Geodetic latitude of observer) [rad].
%          - HA of telescope pole [rad] 
%          - Distance between NCP and equaltorial pole [rad]
% Output : - Tracking error in HA ["/s]
%          - Tracing error in Dec ["/s]
%          - Tracking error in RA ["/s]
%          - HA of telescope [rad]
%          - Dec of telescope [rad]
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [dHt_dt,dDt_dt]=celestial.coo.polar_alignment_drift(1,0,32./RAD,45./RAD,1./RAD)
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

dH_dt = 360.*3600./86164.091; % ["/s]   sidereal rate
dH_dt = dH_dt./(3600.*RAD);

Dt = asin(cos(Beta).*sin(D) + sin(Beta).*cos(D).*cos(Psi - H));

dDt_dt = dH_dt .* (cos(D).*sin(Beta).*sin(Psi-H))./cos(Dt);

Ht = acos(cos(Z)./(cos(Dt).*cos(Phit)) - tan(Dt).*tan(Phit));

%CosHt = cos(Z)./(cos(Dt).*cos(Phit)) - tan(Dt).*tan(Phit);
%Ep = acos((sin(Dt) - cos(Z).*sin(Phit))./(sin(Z).*cos(Phit)));
%SinHt = sin(Ep).*sin(Z)./cos(Dt);
%Ht = atan2(SinHt,CosHt);

[Az,Alt]=celestial.coo.ha2az(H,D,Phi);

if Psi>pi
    % P' is Eastward
    
    
else
    % P' is Westward
    
end

gamma = acos( (cos(Beta) - sin(Phi).*sin(Phit))./(cos(Phi).*cos(Phit)) );



% note that there is an error in the last term of Equation 4 in Markworth
%dHt_dt = dH_dt .* cos(D).*cos(Phi).*sin(H)./(cos(Dt).*cos(Phit).*sin(Ht)) - dDt_dt.*cos(Z).*tan(Dt)./(sin(Ht).*cos(Phit).*cos(Dt)) + ...
%         dDt_dt.*tan(Phit)./(cos(Dt).^2 .* sin(Ht));

% This is incorrect for H=0 and will return NaN

Denom = sin(H).*cos(Dt).*cos(Phi).^2;
dHt_dt = (cos(Phi).*cos(Phit).*cos(D).*cos(Dt).*sin(H).*dH_dt + ...
         sin(Phit).*cos(Phit).*dDt_dt - ...
         cos(Z).*cos(Phit).*sin(Dt).*dDt_dt)./Denom;
     


     
dDt_dt = dDt_dt.*3600.*RAD;
dHt_dt = dHt_dt.*3600.*RAD;
dRt_dt = dHt_dt - dH_dt.*3600.*RAD;
