function Il=fun_lorentzian(Par,Wave)
%------------------------------------------------------------------------------
% fun_lorentzian function                                            AstroSpec
% Description: Calculate a Lorntzian profile spectral line:
%              Cont + D.*Gamma./(pi.*( (X-X0).^2 + Gamma.^2 ))
% Input  : - Vector of parameters [Line Center, Line Width, Line depth].
%          - Vector of wavelength [A] in which to calculate the Lorentzian
%            profile.
%          - Vector (same length as first input argument) or scalar containing
% Output : - Line intensity as a function of input wavelngth [power A^-1]
% Tested : Matlab 7.13
%     By : Eran O. Ofek                        Oct 2012
%    URL : http://www.weizmann.ac.il/home/eofek/matlab/
%------------------------------------------------------------------------------


if (ischar(Par)),
    % get guess for parameters:
    Spec = Wave;
    Par0(1) = mean(Spec(:,2));  % amplitude
    Par0(2) = mean(Spec(:,1));  % central wavelength
    Par0(3) = (max(Spec(:,1)) - min(Spec(:,1)))./10; % std
    
    Il = Par0;
else
    
    Wave0   = Par(1);
    Gamma   = Par(2);
    D       = Par(3);
    
    X  = Wave;
    X0 = Wave0;
 
    Il = D.*Gamma./(pi.*( (X-X0).^2 + Gamma.^2 ));
end
