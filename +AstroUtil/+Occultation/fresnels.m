function FS=fresnels(T)
% Fresnel sine function
% Package: AstroUtil.Occultation
% Description: Return the Fresnel sine function f(T): sin(0.5*pi*T^2)
% Input  : - Parameter.
% Output : - Fresnel sine function.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Apr 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------
FS = sin(0.5.*pi.*T.^2);
