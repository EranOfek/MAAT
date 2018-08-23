function Y=fun_2gauss_conv(Par,WL);

Res = 2;          % Lick red
Res = 2./3.15;    % Lick blue
%Res = 1.25        % KPNO

Res = 0.5;

if (ischar(Par)),
    % get guess for parameters:
    Spec = WL;
    Par0(1) = mean(Spec(:,2));  % amplitude
    Par0(2) = mean(Spec(:,1));  % central wavelength
    Par0(3) = (max(Spec(:,1)) - min(Spec(:,1)))./10; % std
    Par0(4) = Par0(1).*2;
    Par0(5) = Par0(2);
    Par0(6) = Par0(3)./4;
    Par0(7) = Par0(1).*3;
    Par0(8) = Par0(2);
    Par0(9) = Par0(3)./8;
    Y = Par0;
    else
      Par(5) = Par(2);
      Par(8) = Par(2);

   Y = Par(1).*exp(-(WL-Par(2)).^2./(2.*Par(3).^2)) + Par(4).*exp(-(WL-Par(5)).^2./(2.*Par(6).^2)) +Par(7).*exp(-(WL-Par(8)).^2./(2.*Par(9).^2));

   ConvWL = (-10:1:10);
   Res = exp(-ConvWL.^2./(2.*Res.^2));
   Res = Res./trapz(ConvWL,Res);

   Y = conv(Y,Res,'same');
end
