function Y=fun_1gauss(Par,WL);

if (ischar(Par)),
    % get guess for parameters:
    Spec = WL;
    Par0(1) = mean(Spec(:,2));  % amplitude
    Par0(2) = mean(Spec(:,1));  % central wavelength
    Par0(3) = (max(Spec(:,1)) - min(Spec(:,1)))./10; % std
   
    Y = Par0;
else
   Y = Par(1).*exp(-(WL-Par(2)).^2./(2.*Par(3).^2));
end
