function [NewPar,NewParErr,Chi2,Deg,Cov,Resid]=fitpow(X,Y,DelY)
% Fit a power-law function of the form: Y = A * X ^(Gamma)
% Package: Util.fit
% Description: Fit a power-law function of the form: Y = A * X ^(Gamma),
%              to data set.
% Input  : - Vector of independent variable.
%          - Vector of dependent variable.
%          - vector of errors ins dependent variable.
%            if scalar is given then its taken
%            as equal error.
% Output : - vector of parameters [A,Gamma]
%          - vector of errors in parameters [err(A),err(Gamma)]
%          - \chi^2
%          - Degrees of freedom
%          - Covariance matrix
%          - The Y axis residuals vector in log10 space.
%            [calculated error for Chi^2=1,
%            can be calculated from mean(abs(Resid)) ].
% Tested : Matlab 5.0
%     By : Eran O. Ofek                    Nov 1996
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 1
%------------------------------------------------------------------------------
N   = length(X);   % number of observations
Deg =  N - 2;      % degrees of freedom

% building the H matrix
H = [ones(N,1), log(X).*ones(N,1)];

% Linearize the problem:
% NewY = ln(Y) = ln(A) + Gamma * ln(X)

if (length(DelY)==1)
   DelY = DelY.*ones(size(X));
end
 
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

% Transformin Parameter vector to A and Gamma.
NewPar    = [exp(Par(1)); Par(2)];
NewParErr = [NewPar(1).*ParErr(1); ParErr(2)];

%'Number of degree of freedom :'; Deg;
Resid = NewY - H*Par;
Chi2  = sum((Resid./NewYerr).^2);
%'Chi square per deg. of freedom       : '; Chi2/Deg;
%'Chi square error per deg. of freedom : '; sqrt(2/Deg);


Resid = NewY - H*Par;

% plot the data + fit
%errorxy([X,Y,DelY],[1,2,3],'o');
%hold on;
%Np = 100;
%X=[min(X):(max(X)-min(X))./(Np-1):max(X)]';
%NewH = [ones(Np,1), log(X).*ones(Np,1)];
%Yplot=NewH*Par;
%plot(X,exp(Yplot),'r');

