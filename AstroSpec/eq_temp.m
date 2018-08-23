function Teq=eq_temp(R,T,D,A)
%--------------------------------------------------------------------------
% eq_temp function                                               AstroSpec
% Description: Calculate the eqilibrium temperature of a body
%              illuminated by a black-body radiation.
% Input  : - Radiating source radius [cm].
%          - Radiating source temperature [K].
%          - Distance between radiating source and object [cm].
%          - Geometric albedo of object.
% Output : - Equlibrium temperature of object [K].
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Oct 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Teq=eq_temp(696000e5,5700,1.5e13,0.3)
% Reliable: 2
%--------------------------------------------------------------------------

Teq = ((1-A).*R.^2.*T.^4./(4.*D.^2)).^(1./4);
