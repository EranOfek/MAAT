function [AlphaX,AlphaY,A_11,A_22,A_12,Phi]=alpha_sis(Theta,X0,Y0,S,Norm)
% Gravitational deflection for softend spherical isothermal sphere 
% Package: AstroUtil.lensing
% Description: Calculate the gravitational lensing deflection angle for a
%              softend isothermal sphere (SIS).
% Input  : - Image plane position [ThetaX, ThetaY] in which to calculate
%            the deflections
%          - Vector of X0 of point lenses positions.
%          - Vector of Y0 of point lenses positions.
%          - Softend core./core radius
%          - Normaliztaion. In the limit of singular (s=0) and spherical
%            (q=1) model, the normalization is the Einstein radius of the
%            model, and its related to the 1-d velocity dispersion \sigma by:
%            4*pi*(\sigma / c)^2 (Dls/Ds)
% Output : - Deflection of x component (AlphaX). 
%          - Deflection of y component (AlphaY). 
%          - Jacobian matrix A_11 term.
%          - Jacobian matrix A_22 term.
%          - Jacobian matrix A_12 term.
%          - Gravitational potential Phi.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------


X       = Theta(:,1) - X0;
Y       = Theta(:,2) - Y0;
Rad2    = X.^2 + Y.^2;
S2      = S.^2;
%Rad     = sqrt(Rad2); 
SS2R2   = sqrt(S2 + Rad2);
AlphaR_d_Rad  = Norm.*(SS2R2 - S)./Rad2;

AlphaX  = AlphaR_d_Rad .* X;
AlphaY  = AlphaR_d_Rad .* Y; 

if (nargout>2),

   A_11 = (S.*X.^2 + S2.*S + sqrt(S2 + Rad2).*S2+ sqrt(S2 + Rad2).*Y.^2 + S.*Y.^2).*Norm./(S2 + Rad2)./(S + sqrt(S2 + Rad2)).^2;

   A_22 = (S.*X.^2 + S2.*S + sqrt(S2 + Rad2).*S2 + sqrt(S2 + Rad2).*X.^2 + S.*Y.^2).*Norm./(S2 + Rad2)./(S + sqrt(S2 + Rad2)).^2;

   A_12 = -1./sqrt(S2 + Rad2).*Norm.*X.*Y./(S + sqrt(S2 + Rad2)).^2;

end

if (nargout>5),
   Phi     = Norm.*(sqrt(S2 + Rad2) - S.*log( (S + SS2R2 )./(2.*S) ) );
end
