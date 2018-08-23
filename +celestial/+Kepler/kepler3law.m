function Res=kepler3law(Mass,Prop,Val)
% Kepler 3rd law
% Package: celestial.Kepler
% Description: Calculate the velocity, semi-major axis and period of
%              a system using the Kepler third law.
% Input  : - System mass [gr].
%          - Type of given property:
%            'p' - period [s].
%            'a' - semi-major axis [cm].
%            'v' - velocity [cm/s].
%          - Propery value.
% Output : - A structure containing the following fields:
%            .p  - period [s].
%            .a  - semi-major axis [cm].
%            .v  - velocity [cm/s].
% Tested : Matlab R2012a
%     By : Eran O. Ofek                    Nov 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=celestial.Kepler.kepler3law(1.9891e33,'a',1.5e13);
% Reliable: 1
%-----------------------------------------------------------------------------

G=constant.G;

Const = G.*Mass./(4.*pi.^2);

switch lower(Prop)
 case 'p'
    Res.p = Val;
    Res.a = (Const.*Val.^2).^(1./3);
    Res.v = (2.*pi.*G.*Mass./Val).^(1./3);
 case 'a'
    Res.a = Val;
    Res.p = sqrt(Val.^3./Const);
    Res.v = sqrt(G.*Mass./Val);
 case 'v'
    Res.v = Val;
    Res.p = 2.*pi.*G.*Mass./(Val.^3);
    Res.a = G.*Mass./(Val.^2);
 otherwise
    error('Unknown Prop option');
end
