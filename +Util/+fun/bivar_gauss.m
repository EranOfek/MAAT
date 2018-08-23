function F=bivar_gauss(X,Y,Par,Par2)
%------------------------------------------------------------------------------
% bivar_gauss function                                                 General
% Description: Return the value of a normalized bivariate (2-dim)
%              Gaussian: 
%              F = Const + Norm/(2*pi*SigmaX*SigmaY*sqrt(1-Rho^2)) *
%                  exp(-1/(2*(1-Rho^2))*((X-X0)^2/SigmaX^2 +
%                                        (Y-Y0)^2/SigmaY^2 -
%                                        2*Rho*(X-X0)*(Y-Y0)/
%                                        (SigmaX*SigmaY)))
% Input  : - X
%          - Y
%          - Vector of parameters:
%            [X0, Y0, SigmaX, SigmaY, Rho, Norm, Const],
%            or alternatively a structure (Par) containing some of
%            the parameters: 
%            .X0, .Y0, .SigmaX, .SigmaY, .Rho .Norm, .Const
%            Default for Norm is 1.
%            Default for Const is 0.
%          - Optioanl structure (Par2) containing the rest of the
%            Gaussian parameters:
%            (.X0, .Y0, .SigmaX, .SigmaY, .Rho, .Norm, .Const). 
%            Note that Par2 overide Par.
% Output : - Value of bivariate Gaussian at [X,Y].
% See also: chi2_bivar_gauss.m
% Tested : Matlab 7.0
%     By : Eran O. Ofek                       May 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [MatX,MatY] = meshgrid([1:1:10],[1:1:10]);
%          Par.X0 = 5; Par.Y0 = 5;
%          Par2.SigmaX = 2; Par2.SigmaY = 2; Par2.Rho = 0;
%          F=bivar_gauss(MatX,MatY,Par,Par2);
% Reliable: 1
%------------------------------------------------------------------------------
DefConst = 0;
DefNorm  = 1;

if (nargin==3),
   Par2 = [];
end

if (isstruct(Par)==0),
   X0     = Par(:,1);  
   Y0     = Par(:,2);
   SigmaX = Par(:,3);
   SigmaY = Par(:,4);
   Rho    = Par(:,5);
   if (length(Par)==7),
      Norm  = Par(:,6);
      Const = Par(:,7);
   elseif (length(Par)==6),
      Norm  = Par(:,6);
      Const = DefConst
   elseif (length(Par)==5),
      Norm  = DefNorm;
      Const = DefConst;
   else
      error('Illegal parameters');
   end
else
   if (isfield(Par2,'X0')==1),
      X0  = Par2.X0;
   else
      X0  = Par.X0;
   end

   if (isfield(Par2,'Y0')==1),
      Y0  = Par2.Y0;
   else
      Y0  = Par.Y0;
   end

   if (isfield(Par2,'SigmaX')==1),
      SigmaX  = Par2.SigmaX;
   else
      SigmaX  = Par.SigmaX;
   end

   if (isfield(Par2,'SigmaY')==1),
      SigmaY  = Par2.SigmaY;
   else
      SigmaY  = Par.SigmaY;
   end

   if (isfield(Par2,'Rho')==1),
      Rho  = Par2.Rho;
   else
      Rho  = Par.Rho;
   end

   if (isfield(Par2,'Const')==1),
      Const  = Par2.Const;
   elseif (isfield(Par,'Const')==1),
      Const  = Par.Const;
   else   
      Const  = DefConst;
   end

   if (isfield(Par2,'Norm')==1),
      Norm  = Par2.Norm;
   elseif (isfield(Par,'Norm')==1),
      Norm  = Par.Norm;
   else   
      Norm  = DefNorm;
   end
end
   
F = Const + Norm./(2.*pi.*SigmaX.*SigmaY.*sqrt(1-Rho.^2)) .* ...
    exp(-1./(2.*(1-Rho.^2)) .* ...
 	                       ((X-X0).^2./SigmaX.^2 + ...
                                (Y-Y0).^2./SigmaY.^2 - ...
				2.*Rho.*(X-X0).*(Y-Y0)./(SigmaX.*SigmaY)));


