function [Res,Period]=eb_light_curve(Lum,Mass,Radii,OrbEl,LimbFun,LDPars1,LDPars2,RefPars1,RefPars2)
% Eclipsing binary light curve as a function of time.
% Package: AstroUtil.binary
% Description: Calculate eclipsing binary light curve as a function of time.
% Input  : - Stars luminosities, [L1, L2].
%          - Stars mass, [M1, M2]. in solar mass.
%          - Stars radii, [R1, R2]. in au.
%          - Orbital elements : [T, q, e, Om, w, i]
%            Units : days, au, [], radians, radians, radinas respectively.
%          - Limb-darkening function, default is 'limb_darkening'.
%          - Parameters for Primary-star limb-darkening pars.
%            Default is {'Milne',1}.
%          - Parameters for Secondary-star limb-darkening pars.
%            Default is {'Milne',1}.
%          - Parameters for Primary-star reflection-effect pars.
%            Default is [0;0;0].
%          - Parameters for Secondary-star reflection-effect pars.
%            Default is [0;0;0].
% Output : - Light curve, two column matrix [JD, LC].
%            LC is given in relative magnitude.
%          - The period in days, as calculated from the input parameters.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [LC,Period]=AstroUtil.binary.eb_light_curve([100 20],[1 0.6],[0.01 0.01],[0 0.3 0.2 0 0 pi./2],'AstroUtil.binary.limb_darkening');
% Reliable: 2
%------------------------------------------------------------------------------
import AstroUtil.binary.*

Nstep   = 1000;
Def.LimbFun = 'limb_darkening';
Def.LDPars1 = {'Milne',1};
Def.LDPars2 = {'Milne',1};

%NstepLC = 10000;
if (nargin==4)
   LimbFun   = Def.LimbFun;
   LDPars1   = Def.LDPars1;
   LDPars2   = Def.LDPars2;
   RefPars1  = [0;0;0];
   RefPars2  = [0;0;0];
elseif (nargin==5)
   LDPars1   = Def.LDPars1;
   LDPars2   = Def.LDPars2;
   RefPars1  = [0;0;0];
   RefPars2  = [0;0;0];
elseif (nargin==6)
   LDPars2   = Def.LDPars2;
   RefPars1  = [0;0;0];
   RefPars2  = [0;0;0];
elseif (nargin==7)
   RefPars1  = [0;0;0];
   RefPars2  = [0;0;0];
elseif (nargin==8)
   RefPars2  = [0;0;0];
elseif (nargin==9)
   % do nothing
else
   error('Illigal number of input arguments');
end



L1 = Lum(1);
L2 = Lum(2);
R1 = Radii(1);
R2 = Radii(2);


TotL1 = L1.*total_light(R1,LimbFun,Nstep,LDPars1);
TotL2 = L2.*total_light(R2,LimbFun,Nstep,LDPars2);

% read orbital elements
T     = OrbEl(1);
Q     = OrbEl(2);
E     = OrbEl(3);
Omega = OrbEl(4);
OmP   = OrbEl(5);
Inc   = OrbEl(6);
A     = Q./(1-E);

if (Q<=(R1+R2))
   error('Error : q<=R1+R2');
end
if (E>=1)
   error('Error : e>=1');
end

% period relation...
TotM   = Mass(1) + Mass(2);
SY     = 365.256363;           % solar sidereal year
Period = SY.*sqrt(A.^3)./sqrt(TotM);

NstepLC = ceil(100.*(Q./min([R1,R2])));
DelT = Period./NstepLC;
JD = [T:DelT:T+Period].';

%[R,Ni] = kepler(JD,OrbEl,TotM);
K = gauss_grav_const(Mass(2),Mass(1));
[Ni,R] = kepler_elliptic(JD-T, Q, E, K);

XYZ=trueanom2pos(R,Ni,Omega,OmP,Inc);
X = XYZ(:,1);
Y = XYZ(:,2);
Z = XYZ(:,3);

% projected distance
D  = sqrt(X.^2 + Y.^2);
D3 = sqrt(X.^2 + Y.^2 + Z.^2);
% phase angle (observer-primary-secondary)
Ph = atan2(X,Z);

% reflection from L2
Ref2 = binary_reflection_effect(R2,D3,Ph,RefPars2);
% reflection from L1
Ref1 = binary_reflection_effect(R1,D3,pi-Ph,RefPars1);


Izp = find(Z>=0);
Izm = find(Z<0);

% secondary in front
ObsP      = zeros(size(D));
ObsP(Izp) = L1.*obstruction(D(Izp),R1,R2,Nstep,LimbFun,LDPars1);
% primary in front
ObsS      = zeros(size(D));
ObsS(Izm) = L2.*obstruction(D(Izm),R2,R1,Nstep,LimbFun,LDPars2);

LC = (TotL1 - ObsP) + TotL1.*Ref2 + (TotL2 - ObsS) + TotL2.*Ref1;


Res = [JD, -2.5.*log10(LC)];
