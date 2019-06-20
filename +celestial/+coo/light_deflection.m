function P1=light_deflection(P,Q,E)
%------------------------------------------------------------------------------
% light_deflection function                                              ephem
% Description: Calculate the observer-centric direction of a planet,
%              corrected for light deflection in the natural frame.
%              Note: for the stellar case Q=P.
% Input  : - 3 by n matrix (P) of planet observer-centric position vector.
%            3-element column vector per epoch.
%            If more than one column is given, then the light deflection
%            is calculated per column vector. In that case Q and E
%            should contain the same number of columns.
%            Position vector units should be AU.
%          - (Q) Planet helicentric position vector. [au]
%          - (E) Observer helicentric position vector. [au]
% Output : - Observer-centric direction (P1) of the object, corrected
%            for light deflection in the natural frame. [au]
% Tested : Matlab 5.3
%     By : Eran O. Ofek                       May 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: E=[1;0;0]; Q=[0;100;0]; P=[0;99;0];
%          P1=light_deflection(P,Q,E);
%------------------------------------------------------------------------------
MUC2 = 9.8706286e-9;  % mu/c^2 [au] 

AbsP = ones(3,1)*sqrt(sum(P.^2,1));     % P abs. value
AbsQ = ones(3,1)*sqrt(sum(Q.^2,1));     % Q abs. value
AbsE = ones(3,1)*sqrt(sum(E.^2,1));     % E abs. value

% convert to unit vectors
UnitP = P./AbsP;
UnitQ = Q./AbsQ;
UnitE = E./AbsE;

% scalar products
QEsp = ones(3,1)*sum(UnitQ.*UnitE,1);
PQsp = ones(3,1)*sum(UnitP.*UnitQ,1);
EPsp = ones(3,1)*sum(UnitE.*UnitP,1);

P1 = P + (2.*AbsE.*MUC2).*(PQsp.*UnitE - EPsp.*UnitQ)./(1 + QEsp);

P1 = P1./[ones(3,1)*sqrt(sum(P1.*P1))];
