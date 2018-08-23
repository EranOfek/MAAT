function [ThetaE,Theta1,Theta2,Mu1,Mu2,Delay1,Delay2,ThetaC]=pointsource_lens(Mass,D_l,D_s,D_ls,Beta,Type)
% Microleninsg magnification, images positions and time delays (OBSOLETE)
% Package: AstroUtil.microlensing
% Description: Calculate the magnification and images positions for a
%              point source lens, given the lens properties and impact
%              parameter. OBSOLETE, use ps_lens instead.
% Input  : - Deflector mass [solar mass]. 
%          - D_l [pc] or z_l.
%          - D_s [pc] or z_s.
%          - D_ls [pc] or z_ls.
%          - Impact parameter [radians].
%          - 0 (default) if distances are given and 1 for redshift.
% Output : - Einstein radius [radians].
%          - Theta_1 [radians] in respect to the lens.
%          - Theta_2 [radians] in respect to the lens.
%          - Mu_1 
%          - Mu_2
%          - Time delay 1 
%          - Time delay 2 
%          - Theta center of mass of images.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jun 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% [ER,T1,T2,Mu1,Mu2,TD1,TD2,Tcm]=AstroUtil.microlensing.pointsource_lens(1,5000,10000,5000,1./(1000.*RAD.*3600));
%-----------------------------------------------------------------------------
H0     = 70;
OmegaM = 0.3;
OmegaL = 0.7;

if (nargin==5)
   Type = 0;
elseif (nargin==6)
   % do nothing
else
   error('Illegal number of input arguments');
end

G   = constant.G;
C   = constant.c;
Pc  = constant.pc;
SM  = constant.SunM;

switch Type
 case 0
    % do nothing
    Z_l = 0;
 case 1
    Z_l  = D_l;
    D_ls = ad_dist([D_l D_s],H0,OmegaM,OmegaL);
    D_l  = ad_dist(D_l,H0,OmegaM,OmegaL);
    D_s  = ad_dist(D_s,H0,OmegaM,OmegaL);
 otherwise
    error('Unknown Type Option');
end
 
D   = D_ls./(D_l.*D_s);
D   = D./Pc;

ThetaE2 = D.*4.*G.*Mass.*SM./(C.^2);
ThetaE  = sqrt(ThetaE2);                 % radius Einstein [radians]

Theta1 = 0.5.*(Beta + sqrt(Beta.^2 + 4.*ThetaE.^2));
Theta2 = 0.5.*(Beta - sqrt(Beta.^2 + 4.*ThetaE.^2));

U      = Beta./ThetaE;
% magnification
Mu0    = (U.^2 + 2)./(2.*U.*sqrt(U.^2 + 4));
Mu1    = Mu0 + 0.5;
Mu2    = Mu0 - 0.5;


% potential
Phi1    = ThetaE2.*log(abs(Theta1));
Phi2    = ThetaE2.*log(abs(Theta2));


% time delay
Delay1 = (1 + Z_l).*(0.5.*(Theta1 - Beta).^2 - Phi1)./(C.*D);
Delay2 = (1 + Z_l).*(0.5.*(Theta2 - Beta).^2 - Phi2)./(C.*D);

% Theta center of mass
ThetaC = 0.5.*Beta + U.*sqrt(U.^2+4).*sqrt(4.*ThetaE.^2+Beta.^2)./(2.*(U.^2+2));
