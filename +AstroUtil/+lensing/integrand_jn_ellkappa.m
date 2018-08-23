function Val=integrand_jn_ellkappa(U,N,X2,Y2,Q2,Alpha,S2,Norm)
%-----------------------------------------------------------------------------
% integrand_jn_ellkappa function                                        glens
% Description: Calculate the integrand of J_n(x,y), for gravitational
%              lensing softened elliptical mass distribution
%              (See Keeton 2001, Eq. 15).
% Input  : - U, integration parameter
%          - N, function order (e.g., J_n).
%          - X.^2
%          - Y.^2
%          - Q.^2, where Q is the ellipse axis ratio b/a.
%          - Power law slope (Alpha).
%          - S.^2, where S is the core radius.
% Output : - Value of the integrand
% Tested : Matlab 6.5
%     By : Eran O. Ofek                     March 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%-----------------------------------------------------------------------------

% Upsilon^2 (Keeton 2001, Eq. 15)
Ups2 = U.*( X2 + Y2./( 1 - (1-Q2).*U  ) );

% Kappa(upsilon) for elliptical mass
% (Keeton 2001, Eq. 13)
Kappa = 0.5.*Norm.^(2-Alpha)./( (S2 + Ups2).^(1-0.5.*Alpha) );

Val   = Kappa./( ( 1 - (1 - Q2).*U ).^(0.5 + N)  );
