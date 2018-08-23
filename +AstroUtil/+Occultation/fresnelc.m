function FC=fresnelc(T)
% Fresnel cosine function
% Package: AstroUtil.Occultation
% Description: Return the Fresnel cosine function f(T): cos(0.5*pi*T^2)
% Input  : - Parameter.
% Output : - Fresnel cosine function.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Apr 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------
FC = cos(0.5.*pi.*T.^2);
