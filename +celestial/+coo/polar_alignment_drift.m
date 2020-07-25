function [H_dot,D_dot,R_dot,DAz,DAlt]=polar_alignment_drift(H,D,Phi,Psi,Beta)
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
%          - Az needed to move the mount in order to coorect alignement (positive Eastward)
%          - Alt needed to move the mount in order to correct the alignement (positive upward)
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [H_dot,D_dot,R_dot]=celestial.coo.polar_alignment_drift(10./RAD,0,32./RAD,45./RAD,1./RAD)
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


Ht_dot = 360.*3600./86164.091; % ["/s]   sidereal rate

Z      = acos(sin(Phi).*sin(D) + cos(Phi).*cos(D).*cos(H));

Phit   = asin(sin(Phi)*cos(Beta) + cos(Phi).*sin(Beta).*cos(Psi));
Dt     = asin(cos(Beta).*sin(D) + sin(Beta).*cos(D).*cos(Psi - H));
Ht     = acos((cos(Z) - sin(Phit).*sin(Dt))./(cos(Phit).*cos(Dt)));  

SinPsit   = cos(Phi).*sin(Psi)./cos(Phit);
CosPsit   = (sin(Phi) - cos(Beta).*sin(Phit))./(sin(Beta).*cos(Phit));
Psit      = atan2(SinPsit,CosPsit);


[Az,Alt] = celestial.coo.ha2az(H,D,Phi);
% calculate the Az of P'
DAlt = -(Phit - Phi);
DAz  = asin(sin(Beta).*sin(Psi)./cos(Phit));   % Az that need to move the mount in order to get to P (positive Eastward)

% azimuth of P'
AzPt = -DAz;
%if AzPt<0
%    AzPt = 2.*pi + AzPt;
%end


% for large H, H and Ht always have the same sign.
% however, when the star is north of the zenith and its azimuth
% is smaller than DAz the sign may flip
SignFlip = (Az>=0 & Az<(pi./2) & Az<=AzPt) | (Az<=0 & Az>(-pi./2) & Az>=AzPt);
SignFlip = 1-2.*SignFlip;

Ht     = abs(Ht).*sign(H).*SignFlip;   % not accurate near the meridian may flip sign


Hp = acos((cos(Z) - sin(Phi).*sin(D))./(cos(Phi).*cos(D)));

%Ht = -Ht;

D_dot  = -Ht_dot.* sin(Beta).*cos(Dt).*sin(Psit+Ht)./cos(D);

H_dot  = (Ht_dot.*cos(Dt).*cos(Phit).*sin(Ht) + D_dot.*sin(Phit).*cos(D) - D_dot.*sin(D).*cos(Phi).*cos(H))./(cos(Phi).*cos(D).*sin(H));

%H_dot(H>0) = -H_dot(H>0);

%H_dot  = (Ht_dot.*cos(Dt).*cos(Phit).*sin(Ht) + D_dot.*sin(Phit).*cos(D) - D_dot.*sin(D).*cos(Phi).*cos(-H))./(cos(Phi).*cos(D).*sin(-H));

%H_dot  = (Ht_dot.*cos(Dt).*cos(Phit).*sin(Ht) - D_dot.*sin(D).*cos(Phi).*cos(H))./(cos(Phi).*cos(D).*sin(H));

% plot(H,(Ht_dot.*cos(Dt).*cos(Phit).*sin(Ht) )./(cos(Phi).*cos(D).*sin(H)))
% hold on;
% plot(H,(D_dot.*sin(Phit).*cos(D) )./(cos(Phi).*cos(D).*sin(H)))
% plot(H,( -D_dot.*sin(D).*cos(Phi).*cos(H))./(cos(Phi).*cos(D).*sin(H)))


R_dot = -(H_dot - Ht_dot);   % R.A. is measured in opposite direction to H.A.

%D_Dot = -Ht_dot.*sec(D).*cos(pi./2-Beta).*cos(Dt).*sin(Psit);
R_dot = Ht_dot - D_dot.*(tan(D).*sec(pi./2-Beta).*sin(Dt) - sec(D).*tan(pi./2-Beta)).*sec(D).*csc(Psi-H);


if 1==0

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


end