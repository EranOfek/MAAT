function [NewPar,NewParErr,Chi2,Deg,Cov,Resid,ParErrBootstrap]=fitexp(X,Y,DelY);
%------------------------------------------------------------------------------
% fitexp function                                                       FitFun
% Description: Fit an exponential function to data. The function has the
%              form: Y = A * exp(-X./Tau).
% Input  : - Vector of independent variable.
%            Alternatively, a amtrix of [X,Y,DelY].
%            In this case the function expect a single input parameter.
%          - Vector of dependent variable.
%          - Vector of errors in dependent variable.
% Output : - Vector of parameters [A,Tau]
%          - Vector of errors in parameters [err(A),err(Tau)]
%          - \chi^2
%          - Degrees of freedom
%          - Covariance matrix
%          - Vector of residuals
% Tested : Matlab 5.1
%     By : Eran O. Ofek                  November 1996
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------
if (nargin==1),
    DelY = X(:,3);
    Y    = X(:,2);
    X    = X(:,1);
end

N   = length(X);   % number of observations
Deg =  N - 2;      % degrees of freedom

% building the H matrix
H = [ones(N,1), -X.*ones(N,1)];

% Linearize the problem:
% NewY = ln(Y) = ln(A) - X./Tau 
NewY    = log(Y);
NewYerr = DelY./Y;

% The Variance Matrix
V = diag(NewYerr.^2);

% inverse the V Matrix
InvV = inv(V);

% The Covariance Matrix
Cov = inv(H'*InvV*H);

% The parameter vector [ln(A); 1./Tau]
Par    = Cov*H'*InvV*NewY;
ParErr = sqrt(diag(Cov));

% Transformin Parameter vector to A and Tau.
NewPar    = [exp(Par(1)), 1./Par(2)];
NewParErr = [NewPar(1).*ParErr(1), ParErr(2).*NewPar(2).^2];

%'Number of degree of freedom :', Deg
Resid = NewY - H*Par;
Chi2  = sum((Resid./NewYerr).^2);
%'Chi square per deg. of freedom       : ',Chi2/Deg
%'Chi square error per deg. of freedom : ',sqrt(2/Deg)


if (nargout>6),
    [StD]=bootstrap_std([X,Y,DelY],@fitexp,1000);
    ParErrBootstrap = StD;
end

%% plot the data + fit
%%errorxy([X,Y,DelY]);
%hold on;
%Np = 100;
%X=[min(X):(max(X)-min(X))./(Np-1):max(X)]';
%NewH = [ones(Np,1), -X.*ones(Np,1)];
%Yplot=NewH*Par;
%plot(X,exp(Yplot),'r');

