function Lup=luptitude(Flux,Flux0,B)
%--------------------------------------------------------------------------
% luptitude function                                             AstroSpec
% Description: Convert flux to luptitudes (asinh magnitudes).
%              OBSOLETE: Use convert.luptitude instead.
% Input  : - Flux.
%          - Reference flux (Flux0), default is 1.
%            Note that this parameter should equal to 10.^(0.4.*ZP);
%          - Softening parameter (B), default is 1e-10.
% Output : - Luptitude.
% Tested : Matlab 7.13
%     By : Eran O. Ofek                    Jul 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

Def.Flux0 = 1;
Def.B     = 1e-10;
if (nargin==1)
   Flux0 = Def.Flux0;
   B     = Def.B;
elseif (nargin==2)
   B     = Def.B;
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end

Lup = -2.5./log(10).*(asinh((Flux./Flux0)./(2.*B))+log(B));
