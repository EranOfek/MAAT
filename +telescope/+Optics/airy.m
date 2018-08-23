function I_Theta=airy(ApRad,Lambda,Theta,I0)
% Calculate the monochromatic Airy function (circular diffraction)
% Package: telescope.Optics
% Description: Calculate the theoretical diffraction pattern for a perfect
%              circular aperture.
% Input  : - Aperture radius [cm].
%          - Wavelength [cm].
%          - (Vector of) angle for optical axis [radians].
%          - Intensity, at optical axis, default is 1.
% Output : - Intensity at angle theta from the optical axis.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Jan 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%-----------------------------------------------------------------------------

if (nargin==3)
   I0 = 1;
elseif (nargin==4)
   % do nothing
else
   error('Illegal number of input arguments');
end


% The wave number:
K = 2.*pi./Lambda;

KASinT = K.*ApRad.*sin(Theta);

I_Theta = I0.*(2.*besselj(1,KASinT)./KASinT).^2;
