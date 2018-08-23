function [Y,Y1,Y2]=fun_pcyg(Par,WL,Conv
% 1-D simplistic P-Cygni line model
% Package: AstroUtil.spec
% Description: Calculating a P-Cygni profile composed of two Gaussians
%              of the form:
%              Y = AmplitudeEmi*exp( (X-W0)^2/(2*SigmaEmi^2) ) - 
%                  AmplitudeAbs*exp( (X-W0)^2/(2*SigmaAbs^2) ) for X<W0
%                  and only the first term for X>=W0.
%              Optionaly, convolve the the result with a Gaussian.
%              This function can be used by fitting functions like
%              nlinfit_my.m and fit_specline.m
% Input  : - Vector of free parameters of the Gaussian
%            [IntegralEmi, W0, SigmaEmi, HalfIntegralAbs, SigmaAbs, N]
%            Alternatively, if this parameter is 'guess' then the program
%            will attempt to guess (i.e., find initial values for fitting
%            functions) the free parameters, based on the position of
%            the highest value.
%          - Width (in sigma) of convolution kernel.
%            This can be used only if the X position is evenly spaced.
%            Default is 0.
% Output : - The P-Cygni profile function value at X position.
%            Alternatively, if the first input argument is 'guess' then
%            will return the best guess parameters.
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    May 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=[4000:1:9000].'; Y=AstroUtil.spec.fun_pcyg([1 6000 100 0.5 200 2],X,10);
% Reliable: 
%--------------------------------------------------------------------------

Def.Conv = 0;
if (nargin==2)
   Conv = Def.Conv;
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (ischar(Par))
    % get guess for parameters:
    Spec = WL;
    [Max,MaxInd] = max(Spec(:,2));
    [Min,MinInd] = min(Spec(:,2));
    Med          = median(Spec(:,2));
    Par0(1) = Max - Med;
    Par0(2) = Spec(MaxInd,1);
    Par0(3) = abs(Spec(MaxInd,1) - Spec(MinInd,1))./3;
    Par0(4) = Med - Min;
    Par0(5) = Par0(3);
    Par0(6) = 2;
    %Par0(1) = mean(Spec(:,2));  % amplitude
    %Par0(2) = mean(Spec(:,1));  % central wavelength
    %Par0(3) = (max(Spec(:,1)) - min(Spec(:,1)))./10; % std

    Y = Par0;
else
    if (ischar(WL))
       switch lower(WL)
        case 'int'
           % integrate the Gaussian
           
           %Y = sqrt(2.*pi).*Par(1).*Par(3);           
        otherwise
           error('Unknown WL option');
       end
    else
        % construct the P-cygni profile
       Y1 = Par(1).*exp(-(WL-Par(2)).^2./(2.*Par(3).^2))./(sqrt(2.*pi).*Par(3));
       Alpha = Par(6);
       Y2 = exp(-abs(WL-Par(2)).^Alpha./(2.*Par(5).^Alpha))./(sqrt(2.*pi).*Par(5));
       
       %Y2 = Par(4).*exp(-(WL-Par(2)).^2./(2.*Par(5).^2))./(sqrt(2.*pi).*Par(5));
       Y2(WL>Par(2)) = 0;
       Integral = trapz(WL,Y2);
       Y2       = Par(4).*Y2./Integral;
       Y  = Y1 - Y2;
       
       if (Conv>0)
          if (max(diff(WL))~=min(diff(WL)))
             error('If conv>0 data must be evenly spaced');
          end
          Spacing   = mean(diff(WL));
          ConvSigma = Conv./Spacing;
          ConvWL    = (-floor(4.*ConvSigma):1:+ceil(4.*ConvSigma));
          Res       = exp(-ConvWL.^2./(2.*ConvSigma.^2));
          Res       = Res./trapz(ConvWL,Res);

          Y         = conv(Y,Res,'same');
          if (nargout>1)
             Y1        = conv(Y1,Res,'same');
             Y2        = conv(Y2,Res,'same');
          end
      end
   end
end
