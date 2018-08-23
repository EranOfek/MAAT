function PX=confint_probdist(ProbFun,Prob,InterpMethod,Eps);
%------------------------------------------------------------------------------
% confint_probdist function                                          AstroStat
% Description: Calculate two-sided confidence interval of a given numerical
%              probability distribution function.
% Input  : - Two column matrix of the probability distribution function [X,dP]
%          - List of cumulative probabilities in which to report the
%            position X such that P=\int_{-\infty}^{X}{dP(X)dx}.
%            Default is [0.0013, 0.0227, 0.1587, 0.5, 0.8413, 0.9773, 0.9987],
%            which correspinds to -3,-2,-1,0,1,2,3 sigma.
%          - Interpolation method (see interp1.m for options).
%            Default is 'linear'.
%          - Epsilon value (for non monotonic series) default is 1e-12.
% Output : - The X position such that P=\int_{-\infty}^{X}{dP(X)dx}.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: x=[-100:0.1:100]'; p=exp(-x.^2./2);
%          PX=confint_probdist([x,p]);
% Reliable: 2
%------------------------------------------------------------------------------

Def.Prob = [normcdf(0,3,1), normcdf(0,2,1), normcdf(0,1,1), ...
            normcdf(0,0,1), ...
            normcdf(0,-1,1), normcdf(0,-2,1), normcdf(0,-3,1)];
Def.InterpMethod = 'linear';
Def.Eps  = 1e-12;

if (nargin==1),
   Prob = Def.Prob;
   InterpMethod = Def.InterpMethod;
   Eps  = Def.Eps;
elseif (nargin==2),
   InterpMethod = Def.InterpMethod;
   Eps  = Def.Eps;
elseif (nargin==3),
   Eps  = Def.Eps;
elseif (nargin==4),
   % do nothing
else
    error('Illegal number of input arguments');
end


% normalize probability function
X    = ProbFun(:,1);
P    = ProbFun(:,2)./trapz(X,ProbFun(:,2));
P    = (P + Eps);
P    = P./trapz(X,P);
CumP = cumtrapz(X,P);

PX = interp1(CumP,X,Prob,InterpMethod);

