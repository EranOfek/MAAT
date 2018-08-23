function [LC,Period]=plot_eb_lc(Lum,Mass,Radii,OrbEl,LimbFun,LDPars1,LDPars2,RefPars1,RefPars2)
% Plot eclipsing binary light curve as a function of time.
% Package: AstroUtil.binary
% Description: Plot eclipsing binary light curve as a function of time.
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
% Example: AstroUtil.binary.plot_eb_lc([100 20],[1 0.6],[0.01 0.01],[0 0.3 0.2 0 0 pi./2],'AstroUtil.binary.limb_darkening');
% Reliable: 2
%-----------------------------------------------------------------------------
import AstroUtil.binary.*

Def.LimbFun = 'limb_darkening';
Def.LDPars1 = {'Milne', 1};
Def.LDPars2 = {'Milne', 1};

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

   
[LC,Period]=eb_light_curve(Lum,Mass,Radii,OrbEl,LimbFun,LDPars1,LDPars2,RefPars1,RefPars2);

plot(LC(:,1),LC(:,2),'b-');
set(gca,'FontSize',14,'YDir','reverse');
h=xlabel('Days');
set(h,'FontSize',18);
h=ylabel('\Delta{m}');
set(h,'FontSize',18);
