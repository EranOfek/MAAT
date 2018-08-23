1;2cfunction [AlphaX,AlphaY,A_11,A_22,A_12,Phi,Kappa]=calc_alpha(ThetaX,ThetaY,Pars,ModelType);
%-----------------------------------------------------------------------
% calc_alpha function                                             glens
% Description: Call programs for calculating deflection
%                          angles of gravitational lensing.
% Input  : - ThetaX
%          - ThetaY
%          - Model parameters:
%            [X position of center of the mass density (pixels)
%             Y position of center of the mass density (pixels)
%             eccentricity of the potential = sqrt(1-B^2/A^2)
%             position angle of the mass distribution elipse [radians].
%             softend core./core radius
%             Power law
%             Normalization
%             additional parameters]
%          - Model Type:
%            1  - (point) point like mass
%            2  - (spl) softened power law elliptical potential
%            3  - 
%            4  - (SIS) Softend Isothermal Sphere
%            5  - 
%            6  - (kspl) softened power law elliptical density
%            7  -
%            8  -
%            9  - 
% Output : - X Deflection
%          - Y Deflection
%          - Jacobian matrix dAlpha/dTheta _11
%          - Jacobian matrix dAlpha/dTheta _22
%          - Jacobian matrix dAlpha/dTheta _12
% Tested : Matlab 7.0
%     By : Eran O. Ofek      March 2005
%    URL : http://astroclub.tau.ac.il
%-----------------------------------------------------------------------

import AstroUtil.lensing.*

Nmod = length(ModelType);
Nargout = nargout;

if (size(ThetaX,2)==1),
   ThetaVec  = [ThetaX, ThetaY];
   SizeTheta = [NaN 1];
else
   SizeTheta = size(ThetaX);
   %--- Assume ThetaY have the same size ---
   ThetaVec  = [Util.array.mat2vec(ThetaX), Util.array.mat2vec(ThetaY)];
end
Ntheta    = size(ThetaVec,1);

AlphaVecX = zeros(Ntheta,1);
AlphaVecY = zeros(Ntheta,1);
A_11_Vec  = zeros(Ntheta,1);
A_22_Vec  = zeros(Ntheta,1);
A_12_Vec  = zeros(Ntheta,1);

for Imod=1:1:Nmod,

   switch ModelType(Imod)
    case 1
       %--------------------------
       %--- call alpha_point.m ---
       %--------------------------
       switch Nargout
        case 2
           [AlphaX,AlphaY]   = alpha_point(ThetaVec,Pars(Imod,1),Pars(Imod,2),Pars(Imod,7));
        case 5
           [AlphaX,AlphaY,A] = alpha_point(ThetaVec,Pars(Imod,1),Pars(Imod,2),Pars(Imod,7));
           A_11              = A;
           A_22              = A;
           A_12              = zeros(size(A));
        otherwise
           error('Illegal number of ourput arguments');
       end

    case 2
       %------------------------
       %--- call alpha_spl.m ---
       %------------------------
       switch Nargout
        case 2
           [AlphaX,AlphaY]                          = alpha_spl(ThetaVec,Pars(Imod,1),Pars(Imod,2),Pars(Imod,3),Pars(Imod,4),Pars(Imod,5),Pars(Imod,6),Pars(Imod,7));
        case 5
           [AlphaX,AlphaY,A_11,A_22,A_12]           = alpha_spl(ThetaVec,Pars(Imod,1),Pars(Imod,2),Pars(Imod,3),Pars(Imod,4),Pars(Imod,5),Pars(Imod,6),Pars(Imod,7));
        case 6
           [AlphaX,AlphaY,A_11,A_22,A_12,Phi]       = alpha_spl(ThetaVec,Pars(Imod,1),Pars(Imod,2),Pars(Imod,3),Pars(Imod,4),Pars(Imod,5),Pars(Imod,6),Pars(Imod,7));
        case 7
           [AlphaX,AlphaY,A_11,A_22,A_12,Phi,Kappa] = alpha_spl(Theta,Pars(Imod,1),Pars(Imod,2),Pars(Imod,3),Pars(Imod,4),Pars(Imod,5),Pars(Imod,6),Pars(Imod,7));
        otherwise
           error('Illegal number of ourput arguments');
       end

    case 4
       %------------------------
       %--- call alpha_sis.m ---
       %------------------------
       switch Nargout
        case 2
           [AlphaX,AlphaY]                    = alpha_sis(ThetaVec,Pars(Imod,1),Pars(Imod,2),Pars(Imod,5),Pars(Imod,7));
        case 5
           [AlphaX,AlphaY,A_11,A_22,A_12]     = alpha_sis(ThetaVec,Pars(Imod,1),Pars(Imod,2),Pars(Imod,5),Pars(Imod,7));
        case 6
           [AlphaX,AlphaY,A_11,A_22,A_12,Phi] = alpha_sis(ThetaVec,Pars(Imod,1),Pars(Imod,2),Pars(Imod,5),Pars(Imod,7));

        otherwise
           error('Illegal number of ourput arguments');
       end

    case 6
       %------------------------------
       %--- call alpha_kspl_fast.m ---
       %------------------------------
       switch Nargout
        case 2
           [AlphaX,AlphaY]                = alpha_kspl_fast(ThetaVec,Pars(Imod,1),Pars(Imod,2),Pars(Imod,3),Pars(Imod,4),Pars(Imod,5),Pars(Imod,6),Pars(Imod,7));
        case 5
           [AlphaX,AlphaY,A_11,A_22,A_12] = alpha_kspl_fast(ThetaVec,Pars(Imod,1),Pars(Imod,2),Pars(Imod,3),Pars(Imod,4),Pars(Imod,5),Pars(Imod,6),Pars(Imod,7));
        otherwise
           error('Illegal number of ourput arguments');
       end
    otherwise
       error('Unknown ModelType');
   end
   AlphaVecX = AlphaVecX + AlphaX;
   AlphaVecY = AlphaVecY + AlphaY;


   if (nargout>2),
      A_11_Vec  = A_11_Vec  + A_11;
      A_22_Vec  = A_22_Vec  + A_22;
      A_12_Vec  = A_12_Vec  + A_12;
   end
end

if (size(ThetaX,2)==1),
   AlphaX = AlphaVecX;
   AlphaY = AlphaVecY;
   if (nargout>2),
     A_11   = A_11_Vec;
     A_22   = A_22_Vec;
     A_12   = A_12_Vec;
   end
else
   AlphaX = vec2mat(AlphaVecX,SizeTheta(2)).';
   AlphaY = vec2mat(AlphaVecY,SizeTheta(2)).';
   if (nargout>2),
      A_11   = vec2mat(A_11_Vec,SizeTheta(2)).';
      A_22   = vec2mat(A_22_Vec,SizeTheta(2)).';
      A_12   = vec2mat(A_12_Vec,SizeTheta(2)).';
   end
end


