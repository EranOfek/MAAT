function [AlphaX,AlphaY,A_11,A_22,A_12]=alpha_kspl_fast(Theta,X0,Y0,Eccen,PA,S,Alpha,Norm)
% Deflection for softened power law elliptical density
% Package: AstroUtil.lensing
% Description: Calculate deflection for softened power law elliptical
%              density of the form:
%              kappa = 0.5*b^(2-alpha)/((s^2 + r^2)^(1-0.5*alpha))
% Input  : - position in which to calculate deflection and
%            magnification [Theta_x,Theta_y], in pixels.
%          - X position of center of the mass density (pixels)
%          - Y position of center of the mass density (pixels)
%          - eccentricity of the potential = sqrt(1-B^2/A^2)
%          - position angle of the mass distribution elipse [radians].
%          - softend core.
%          - Power law
%          - Normalization
% Output : - X Deflection
%          - Y Deflection
%          - Jacobian matrix dAlpha/dTheta _11
%          - Jacobian matrix dAlpha/dTheta _22
%          - Jacobian matrix dAlpha/dTheta _12
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Mar 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference : Keeton 2002 p.7
% Example: [AlphaX,AlphaY] = alpha_gnfw([ThetaX, ThetaY],...
%                                ModelPars(Icomp,ColX0),...
%                                ModelPars(Icomp,ColY0),...
%                                ModelPars(Icomp,ColE),...
%                                ModelPars(Icomp,ColPA),...
%                                ModelPars(Icomp,ColS),...
%                                ModelPars(Icomp,ColGamma),...
%                                ModelPars(Icomp,ColNorm));
% Needed : compile alpha_kspl_fast_mp.f  alpha_kspl_fast_sp.f
%          mex alpha_kspl_fast_mp.f
%          mex alpha_kspl_fast_sp.f
%------------------------------------------------------------------------------

import AstroUtil.lensing.*

%  [AlphaX,AlphaY]=alpha_kspl_fast([1 1],-1,-5,0.7,1,1,0.8,1)

X          = Theta(:,1) - X0;
Y          = Theta(:,2) - Y0;

% turn the coordinates so that its major axis is on x axis.
X_1        = X.*cos(-PA) - Y.*sin(-PA);
Y_1        = X.*sin(-PA) + Y.*cos(-PA);

S2         = S.^2;
Q2         = 1 - Eccen.^2;
AxisRatio  = sqrt(Q2);


Gamma      = 1 - 0.5.*Alpha;
N          = size(Theta,1);
if (length(AxisRatio)==1 & length(S2)==1 & length(Alpha)==1),
   Pars = [Gamma,AxisRatio,S2];
   AlphaMag = alpha_kspl_fast_sp([X_1, Y_1].',Pars);
else
   Pars = [Gamma,AxisRatio,S2];
   AlphaMag = alpha_kspl_fast_mp([X_1, Y_1].',Pars.');
end

ColX   = 1;
ColY   = 2;
AlphaX = Norm.*(AlphaMag(ColX,:).'.*cos(PA) - AlphaMag(ColY,:).'.*sin(PA));
AlphaY = Norm.*(AlphaMag(ColX,:).'.*sin(PA) + AlphaMag(ColY,:).'.*cos(PA));

if (nargout>2),
   A_11_r   = Norm.*AlphaMag(3,:);
   A_22_r   = Norm.*AlphaMag(6,:);
   A_12_r   = Norm.*AlphaMag(4,:);

   %--- rotate TAT^t
   CosP2 = cos(PA).^2;
   SinP2 = sin(PA).^2;
   Cos2P = cos(2.*PA);
   Sin2P = sin(2.*PA);
   CosSin= cos(PA).*sin(PA);
   A_11  = A_11_r.*CosP2  - A_12_r.*Sin2P + A_22_r.*SinP2;
   A_22  = A_11_r.*SinP2  + A_12_r.*Sin2P + A_22_r.*CosP2;
   A_12  = A_11_r.*CosSin + A_12_r.*Cos2P - A_22_r.*CosSin;

   A_11 = A_11.';
   A_22 = A_22.';
   A_12 = A_12.';

end


