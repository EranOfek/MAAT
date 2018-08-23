function [AlphaX,AlphaY,A_11,A_22,A_12,Phi,Kappa]=alpha_spl(Theta,X0,Y0,e,PA,S,Gamma,Norm)
% Gravitational deflection of softened power law elliptical potential
% Package: AsttroUtil.lensing
% Description: Calculate gravitational lensing deflection and magnification
%              tensor for softened power law elliptical potential of the
%              form: phi = b(s^2 + x^2 + y^2/q^2)^(alpha/2) - b/s^alpha
% Input  : - position in which to calculate deflection and
%            magnification [Theta_x,Theta_y], in pixels.
%          - X position of center of the mass density (pixels)
%          - Y position of center of the mass density (pixels)
%          - eccentricity of the potential = 1-sqrt(1-B^2/A^2)
%          - position angle of the mass distribution elipse [radians].
%          - softend core.
%          - Power law
%          - Normalization
% Output : - X Deflection
%          - Y Deflection
%          - Mapping matrix term A_11
%          - Mapping matrix term A_12
%          - Mapping matrix term A_12
%          - The potential
%          - Kappa, the two dimensional density
% Tested : Matlab 6.5
%     By : Eran O. Ofek & Keren Sharon     Mar 2005
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
% Reliable: 2
%------------------------------------------------------------------------------


InvQ = 1./sqrt(1- e.^2);

% turn the coordinates so that its major axis is on x axis.
X1 =  (Theta(:,1)-X0).*cos(-PA) - (Theta(:,2)-Y0).*sin(-PA);
Y1 = ((Theta(:,1)-X0).*sin(-PA) + (Theta(:,2)-Y0).*cos(-PA)).*InvQ; % stretch elipse to a circle

%--- Calculate Alpha ---
Temp     = Norm.*Gamma.*(S.^2 + X1.^2 + Y1.^2).^(0.5.*Gamma - 1);
AlphaX_r = Temp.*X1;
AlphaY_r = Temp.*Y1;


% turn back to the original position angle
AlphaX = AlphaX_r.*cos(PA) - AlphaY_r.*sin(PA);
AlphaY = AlphaX_r.*sin(PA) + AlphaY_r.*cos(PA);

if (nargout>2),
   TempXY = S.^2 + X1.^2 + Y1.^2;
   Temp0 = Norm.*TempXY.^(0.5.*Gamma - 2).*Gamma;
   Temp1 = Temp0.*Gamma;
   Temp2 =  Norm.*TempXY.^(0.5.*Gamma - 1).*Gamma;
   %A11
   A_11_r = Temp1.*X1.^2 + Temp2 - 2.*Temp0.*X1.^2;
   %A22
   A_22_r = Temp1.*Y1.^2 + Temp2 - 2.*Temp0.*Y1.^2;
   %A12
   A_12_r = Temp1.*Y1.*X1 - 2.*Temp0.*X1.*Y1;

   %--- rotate TAT^t

   CosP2 = cos(PA).^2;
   SinP2 = sin(PA).^2;
   Cos2P = cos(2.*PA);
   Sin2P = sin(2.*PA);
   CosSin= cos(PA).*sin(PA);
   A_11  = A_11_r.*CosP2  - A_12_r.*Sin2P + A_22_r.*SinP2;
   A_22  = A_11_r.*SinP2  + A_12_r.*Sin2P + A_22_r.*CosP2;
   A_12  = A_11_r.*CosSin + A_12_r.*Cos2P - A_22_r.*CosSin;


   if (nargout>5),
      %--- Calculate the potential ---
      Phi = Norm.*TempXY.^(0.5.*Gamma) - Norm/(S.^Gamma);

      if (nargout>6),
         %--- Calculate Kappa ---

         %--- not normalized ---
         error('Not normalized');
         Kappa = 1./Gamma .* (TempXY).^(0.5.*Gamma-1) .* ( (Gamma - 2).*(X1.^2 + (Y1.*InvQ.^2).^2)./TempXY + (1 + InvQ.^2) );
      end
   end

end




