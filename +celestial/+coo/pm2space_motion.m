function [SpaceMotion,SpaceVec]=pm2space_motion(RA,Dec,PM_RA,PM_Dec,Parallax,RadVel)
%--------------------------------------------------------------------------
% pm2space_motion function                                           ephem
% Description: Convert proper motion, radial velocity and parralax to
%              space motion vector in the equatorial system.
% Input  : - [RA] in radians.
%          - [Dec] in radians.
%          - PM in RA [mas/yr]
%          - PM in Dec [mas/yr]
%          - Parallax [mas]
%          - Radial velocity [km/s]
% Output : - Space motion vector [X Y Z] in au/day, in the equatorial
%            system.
%          - Space position vector in au, in the equatorial system. 
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: Seidelmann 1992 p.121
% Reliable: 2
%--------------------------------------------------------------------------
RAD  = 180./pi;
S    = 2.*pi./(360.*3600.*36525);          % convert from "/cy to radians/day
AU   = constant.au('SI')./1000;      % [km]
K    = 86400./AU;                          % convert from km/s to au/day

% conver PM_RA from mas/yr to second of time per cy
PM_RA  = PM_RA .* 100./1000 ./15 ./cos(Dec);

% conver PM_Dec from mas/yr to arcsecond  per cy
PM_Dec = PM_Dec.* 100./1000;

% make sure parallax is not zero
Parallax(abs(Parallax)<1e-4) = 1e-4;   % mas

% Convert Parallax from mas to radians
Parallax = Parallax ./(1000.*3600.*RAD);

N = max([length(RA);length(Dec);length(PM_RA);length(PM_Dec);length(Parallax);length(RadVel)]);
RA      = RA.*ones(N,1);
Dec     = Dec.*ones(N,1);
PM_RA   = PM_RA.*ones(N,1);
PM_Dec  = PM_Dec.*ones(N,1);
Parallax= Parallax.*ones(N,1);
RadVel  = RadVel.*ones(N,1);

R = 1./sin(Parallax);


SpaceVec = celestial.coo.cosined([RA,Dec]);
SpaceVec = bsxfun(@times,SpaceVec,R);

SpaceMotion = zeros(N,3);
for I=1:1:N

   Mat = [-cos(Dec(I)).*sin(RA(I)), -sin(Dec(I)).*cos(RA(I)), cos(Dec(I)).*cos(RA(I)); ...
           cos(Dec(I)).*cos(RA(I)), -sin(Dec(I)).*sin(RA(I)), cos(Dec(I)).*sin(RA(I)); ...
	   0                      ,  cos(Dec(I))            , sin(Dec(I))];

   Vec = [15.*S.*R(I).*PM_RA(I); S.*R(I).*PM_Dec(I); K.*RadVel(I)];

   SpaceMotion(I,:) = (Mat*Vec).';
end


