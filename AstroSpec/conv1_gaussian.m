function [Yconv,X]=conv1_gaussian(Y,Sigma,varargin)
%--------------------------------------------------------------------------
% conv1_gaussian function                                        AstroSpec
% Description: Convolve a 1-D signal with a Gaussian kernel.
% Input  : - A column vector of evenly spaced signal (Y).
%            In this case the program assumes that the the series is
%            spaced with unit steps.
%            Alternatively, this can be a two column matrix of [X,Y],
%          - Sigma of the Gaussian to convolve, in units of of the
%            space between elements in Y. If X is not provided, then
%            the function assumes that diff(X)=1, while if X is provided,
%            then Sigma is in the units specified in X.
%            Default is 1.
%          * Arbitrary number of pairs of input arguments ...,key,val,...
%            The following keywords are available:
%            'Bounds'   - Lower and upper range of the Gaussian kernel in
%                         units of of the space between elements in Y.
%                         Default is [floor(-3.*Sigma), ceil(3.*Sigma)].
%                         If empty matrix then use default.
%            'Shape'    - The resultant convolution size
%                         {'full'|'same'|'valid'}.
%                         See conv.m for details. Default is 'same'.
%                         If empty matrix use default.
% Output : - The convolved signal.
%          - X vector.
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Jan 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Y = zeros(100,1); Y(50)=10; YConv=conv1_gaussian(Y,10);
%          X = [1000:1:10000]'; Y=zeros(size(X)); Y(50)=2; Y(7000)=3;
%          [YConv,X]=conv1_gaussian([X,Y],10);
%          YConv=conv1_gaussian([X,Y],10,'Bounds',[-40 40],'Shape','same');
% Reliable: 2
%--------------------------------------------------------------------------

Def.Sigma  = 1;
if (nargin==1),
    Sigma  = Def.Sigma;
else
    % do nothing
end

DefV.Bounds = [];
DefV.Shape  = 'same';
InPar = set_varargin_keyval(DefV,'y','def',varargin{:});

if (isempty(InPar.Bounds)),
    InPar.Bounds  = [-3.*Sigma, +3.*Sigma];
end


if (size(Y,2)==1),
    N = length(Y);
    X = (1:1:N)';
else
    N = size(Y,1);
    X = Y(:,1);
    Y = Y(:,2);
end

if (range(diff(X))>0),
    error('Series required to be evenly spaced');
end

StepSize = min(diff(X));

Sigma          = Sigma./StepSize;
InPar.Bounds   = InPar.Bounds./StepSize;

InPar.Space = 'linear';  % log is TBD
switch lower(InPar.Space)
    case 'linear'
        ConvW = (floor(min(InPar.Bounds)):1:ceil(max(InPar.Bounds)));
        Gauss = exp(-ConvW.^2./(2.*Sigma.^2));
        Gauss = Gauss./trapz(ConvW, Gauss);

        Yconv = conv(Y,Gauss,InPar.Shape);
    case 'log'
        %InPar.InterpMethod = 'linear';
        
        %ConvW = (floor(min(InPar.Bounds)):1:ceil(max(InPar.Bounds)));
        %Gauss = exp(-ConvW.^2./(2.*Sigma.^2));

        %LogX    = log10(X);
        %DiffX   = min(diff(X)./X(1:end-1))./10;
        
        %LogXvec = (min(LogX):DiffX:max(LogX)).';
        %Xvec    = 10.^LogXvec;
        
        %Ylog    = interp1(X,Y,Xvec,InPar.InterpMethod);
        %Yconv = conv(Ylog,Gauss,InPar.Shape);
        %X = Xvec;
    otherwise
        error('Unknown Space option');
end

