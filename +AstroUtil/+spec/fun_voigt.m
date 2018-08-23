function Il=fun_voigt(Par,Wave)
% 1-D Voight profile
% Package: AstroUtil.spec
% Description: Calculate a Voigt line profile (convolution of a Gaussian
%                and a Lorntzian).
% Input  : - Vector of parameters [Line center, Line width, Line depth,
%            Gaussian sigma].
%          - Vector of wavelength [A] in which to calculate the Voigt
%            profile.
% Output : - Line intensity as a function of input wavelngth [power A^-1]
% Tested : Matlab 7.13
%     By : Eran O. Ofek                        Oct 2012
%    URL : http://www.weizmann.ac.il/home/eofek/matlab/
%------------------------------------------------------------------------------


if (ischar(Par))
    % get guess for parameters:
    Spec = Wave;
    Par0(1) = mean(Spec(:,1));  % central wavelength
    Par0(2) = (max(Spec(:,1)) - min(Spec(:,1)))./10; % std
    Par0(3) = mean(Spec(:,2));  % amplitude
    Par0(4) = (max(Spec(:,1)) - min(Spec(:,1)))./20; % std
    
    Il = Par0;
else
    Wave0 = Par(1);
    Gamma = Par(2);
    D     = Par(3);
    Sigma = Par(4);
    
    Gauss = exp(-0.5.*((Wave-Wave0)./Sigma).^2)./(Sigma.*sqrt(2.*pi));

    Il = AstroUtil.spec.fun_lorentzian([Wave0,Gamma,D],Wave);
    Il = conv(Il,Gauss,'same');
end
