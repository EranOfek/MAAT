function Res=fit_gauss1da(X,Y,Err)
% Fit a 1-D Gaussian without background using LSQ.
% Package: Util.fit
% Description: Fit a 1-D Gaussian without background using LSQ.
%              The Gaussian form is:
%              Y=A/(Sigma*sqrt(2*pi))*exp(-(X-X0).^2./(2.*Sigma^2))
% Input  : - X
%          - Y
%          - Error in Y. Default is 1.
% Output : - Structure containing the results.
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Apr 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=[4000:1:9000].'; S=5;
%          Y=10./(sqrt(2.*pi).*S).*exp(-(X-6000).^2./(2.*S.^2));
%          Res=Util.fit.fit_gauss1da(X,Y)
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==2),
    Err = 1;
end


% ln y = ln A - ln(Sigma*sqrt(2*pi))-1/(2.*Sigma^2)*(X^2 - 2.*X.*X0 +X0.^2)
N = numel(X);
Err = Err.*ones(N,1);

FlagGood = Y>0.00001;
X   = X(FlagGood);
Y   = Y(FlagGood);
Err = Err(FlagGood); 
N = numel(X);

% 1st par: ln(A) - ln(Sigma.*sqrt(2.*pi)) - X0.^2./(2.*Sigma.^2)
% 2nd par: +2.*X0./(2.*Sigma.^2)
% 3rd par: -1./(2.*Sigma.^2)

H = [ones(N,1),  X, X.^2];
[Par,ParErr] = lscov(H,log(Y),1./((Err./Y).^2));
Res.X0    = -Par(2)./(2.*Par(3));
Res.Sigma = (sqrt(2).*(-1./Par(3)).^(1./2))./2;
Res.A     = pi.^(1./2).*exp(-Par(2).^2./(4.*Par(3)) + Par(1)).*(-1./Par(3)).^(1./2);
Res.Resid   = Y - Res.A./(Res.Sigma.*sqrt(2.*pi)) .* exp(-(X-Res.X0).^2./(2.*Res.Sigma.^2) );
Res.Chi2    = sum((Res.Resid./Err).^2);
Res.Nobs    = length(X);
Res.Ndof    = Res.Nobs - 3;

