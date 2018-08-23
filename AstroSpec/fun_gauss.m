function Y=fun_gauss(Par,WL,Conv)
%------------------------------------------------------------------------------
% fun_gauss function                                                 AstroSpec
% Description: Calculating a Gaussian function of the form:
%              Y = Amplitude*exp( (X-W0)^2/(2*Sigma^2) )
%              Optionaly, convolve the the result with a Gaussian.
%              This function can be used by fitting functions like
%              nlinfit_my.m and fit_specline.m
% Input  : - Vector of free parameters of the Gaussian
%            [Integral, W0, Sigma]
%            Alternatively, if this parameter is 'guess' then the program
%            will attempt to guess (i.e., find initial values for fitting
%            functions) the free parameters, based on the position of
%            the highest value.
%          - X positions at which to calculate the Gaussian.
%            If the first argument is 'guess' then this parameter should
%            be [X,Y], and the program will attempt to guess the free
%            parameters.
%            Alternatively, this parameter can be 'int' and in this case
%            the program will calculate the integral of the Gaussian.
%          - Width (in sigma) of convolution kernel.
%            This can be used only if the X position is evenly spaced.
%            Default is 0.
% Output : - The Gaussian value at X position.
%            Alternatively, if the first input argument is 'guess' then
%            will return the best guess parameters.
% Tested : Matlab 2011b
%     By : Eran O. Ofek        Apr 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=[4000:1:9000].'; Y=fun_gauss([1 6000 100],X);
%          Y=fun_gauss([1 6000 100],X,100);
%          GyessPar=fun_gauss('guess',[X Y]); % find initial parameters
%          Int = fun_gauss([1 6000 100],'int')
% Reliable: 2
%------------------------------------------------------------------------------

Def.Conv = 0;
if (nargin==2),
   Conv = Def.Conv;
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end


%Res = 2;          % Lick red
%Res = 2./3.15;    % Lick blue
%Res = 1.25        % KPNO
%Res = 2.06;       % Gemini blue
%Res = 3.09;       % Gemini red

if (ischar(Par)),
    % get guess for parameters:
    Spec = WL;
    [Max,MaxInd] = max(Spec(:,2));
    [Min,MinInd] = min(Spec(:,2));
    Par0(1) = Max-Min;
    Par0(2) = Spec(MaxInd,1);
    Par0(3) = abs(Spec(MaxInd,1) - Spec(MinInd,1))./3;

    %Par0(1) = mean(Spec(:,2));  % amplitude
    %Par0(2) = mean(Spec(:,1));  % central wavelength
    %Par0(3) = (max(Spec(:,1)) - min(Spec(:,1)))./10; % std

    Y = Par0;
else
    if (ischar(WL)),
       switch lower(WL)
        case 'int'
           % integrate the Gaussian
           %Y = sqrt(2.*pi).*Par(1).*Par(3); 
           Y = Par(1);
        otherwise
           error('Unknown WL option');
       end
    else
       Y = Par(1).*exp(-(WL-Par(2)).^2./(2.*Par(3).^2))./(sqrt(2.*pi).*Par(3));

       if (Conv>0),
          if (max(diff(WL))~=min(diff(WL))),
 	     error('If conv>0 data must be evenly spaced');
          end
          Spacing   = mean(diff(WL));
          ConvSigma = Conv./Spacing;
          ConvWL    = (-floor(4.*ConvSigma):1:+ceil(4.*ConvSigma));
          Res       = exp(-ConvWL.^2./(2.*ConvSigma.^2));
          Res       = Res./trapz(ConvWL,Res);

          Y         = conv(Y,Res,'same');
      end
   end
end
