function Res=fit_gauss1d(X,Y,Err)
% Fit a 1-D Gaussian without background using fminsearch.
% Package: Util.fit
% Description: Fit a 1-D Gaussian without background using fminsearch.
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
%          Res=fit_gauss1d(X,Y)
% Reliable: 2
%--------------------------------------------------------------------------
import Util.fit.*

if (nargin==2)
    Err = 1;
end

% make sure that the data is double as fminsearch requires double 
X   = double(X);
Y   = double(Y);
Err = double(Err);

%A = Par(1)
%Sigma = Par(2)
%X0    = Par(3);

Fun = @(Par,X,Y,Err) sum((Y - Par(1)./(Par(2).*sqrt(2.*pi)) .* exp(-(X-Par(3)).^2./(2.*Par(2).^2) ) ).^2  ./(Err.^2));
[~,MaxI]   = max(Y);
GuessX0    = X(MaxI);
GuessA     = trapz(X,Y);
GuessSigma = sqrt(sum(Y.*(X-GuessX0).^2)./sum(Y));

Options = optimset(@fminsearch);
Options.MaxFunEvals = 1e5;
Options.MaxIter     = 1e3;
[BestPar,Chi2]=fminsearch_my({Fun,X,Y,Err},[GuessA GuessSigma GuessX0],Options);

Res.A       = BestPar(1);
Res.Sigma   = BestPar(2);
Res.X0      = BestPar(3);
Res.Resid   = Y - Res.A./(Res.Sigma.*sqrt(2.*pi)) .* exp(-(X-Res.X0).^2./(2.*Res.Sigma.^2) );
Res.Chi2    = Chi2;
Res.Nobs    = length(X);
Res.Ndof    = Res.Nobs - length(BestPar);


%plot(X,Y,'k-');
%hold on
%plot(X,Y+Res.Resid,'r-')



