function NewY=conv_gauss1d(X,Y,Sigma,HalfSize,Space,RefX);
%----------------------------------------------------------------------------
% conv_gauss1d function                                               ImSpec
% Description: Convolve a 1-dimension vector with a gaussian.
% Input  : - Equally spaced X of vector to convolve.
%          - Y of vector to convolve.
%          - Sigma of Gaussian to convolve (FWHM = 2.35 * Sigma).
%          - Half window size of convoltion kernel. Outside the window,
%            the effective kerenel is 0.
%          - Space type to convolve in:
%            {'linear' | 'log'}, default is 'linear'.
%          - Reference X for which the Sigma refers to.
% Output : - Convolved Y.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                      July 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%----------------------------------------------------------------------------
OverSampling   = 10;   % oversampling for log space
InterpMethod   = 'linear';

DefSpace = 'linear';
DefRefX  = NaN;
if (nargin==4),
   Space = DefSpace;
   RefX  = DefRefX;
else
   % do nothing
end

switch lower(Space)
 case 'log'
    error('log Space convoltion not working yet!');
    OrigX        = X;
    OrigY        = Y;
    OrigHalfSize = HalfSize;
    OrigSigma    = Sigma;

    X     = log(X);
    Nx    = length(X);

    HSS   = ceil(HalfSize./Sigma);
    Sigma = Sigma./RefX;    % derivative of log(x)

    HalfSize = Sigma.*HSS;

    % resample in log space
    MinX = min(X);
    MaxX = max(X);
    DX   = range(X)./(Nx.*OverSampling);
    X    = [MinX:DX:MaxX].';
    Y    = interp1(OrigX,OrigY,exp(X),InterpMethod);
    GX = [-HalfSize:DX:HalfSize].';
    GY = (1./(Sigma.*sqrt(2.*pi))).*exp(-GX.^2./(2.*Sigma.^2));

    NewY = conv(Y,GY);
    Nny  = length(NewY);

    NewY = NewY([round(HalfSize./DX)+1:1:Nny-round(HalfSize./DX)]);


%    X    = exp(X);
trapz(X,NewY)
    NewY = interp1(X,NewY,OrigX,InterpMethod,'extrap');


 case {'linear','lin'}

    DX = min(diff(X));

    GX = [-HalfSize:DX:HalfSize].';
    GY = (1./(Sigma.*sqrt(2.*pi))).*exp(-GX.^2./(2.*Sigma.^2));

    NewY = conv(Y,GY);
    Nny  = length(NewY);

    NewY = NewY([round(HalfSize./DX)+1:1:Nny-round(HalfSize./DX)]);

 otherwise
    % do nothing
end





