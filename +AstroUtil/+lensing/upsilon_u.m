function Val=upsilon_u(U,X,Y,Q)
%-----------------------------------------------------------------------------
% upsilon_u function                                                    glens
% Description: Calculate the upsilon(u) function (Eq. 15 in Keeton 2001).
% Input  : - u
%          - x
%          - y
%          - q, the prjected axis ratio = b/a
% Output : - Function value
% Tested : Matlab 6.5
%     By : Eran O. Ofek                     March 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%-----------------------------------------------------------------------------

Val = U.*( X.^2 + Y.^2./( 1 - (1-Q.^2).*U  ) );
Val = sqrt(Val);
