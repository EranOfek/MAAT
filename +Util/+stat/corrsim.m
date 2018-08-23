function [Corr,Prob,SimCorr,N]=corrsim(X,Y,Nsim,PermType,varargin)
% Correlation between two vectors and confidence region using bootstrap
% Package: Util.stat
% Description: Calculate the correlation between two vectors and use the
%              bootstrap method to estimate the probability to get
%              a correlation larger than the observed correlation.
%              The function ignores NaN values.
% Input  : - X vector.
%          - Y vector.
%          - Number of bootstrap simulations. Default is 1000.
%          - Which vector to permute in bootstrap test {'x'|'y'}.
%            Default is 'y'.
%          * Arbitrary number of input argument to pass to the
%            built in function corr.m
%            'type' - 'Pearson' (default) to compute Pearson's linear
%                     correlation coefficient.
%                     'Kendall' to compute Kendall's tau.
%                     'Spearman' to compute Spearman's rho.
% Output : - The Correlation.
%          - Probability to get larger correlation.
%          - Vector of simulated correlations.
%          - Number of entries used to calculate correlation.
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Mar 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Corr,Prob,SimCorr,N]=corrsim(rand(100,1),rand(100,1));
% Reliable: 1
%------------------------------------------------------------------------------


Def.Nsim     = 1000;
Def.PermType = 'Y';

if (nargin==2)
   Nsim  = Def.Nsim;
   PermType = Def.PermType;
elseif (nargin==3)
   PermType = Def.PermType;
else
   % do nothing
end

if (isempty(Nsim))
   Nsim  = Def.Nsim;
end
if (isempty(PermType))
   PermType = Def.PermType;
end

if (size(X,2)>size(X,1))
   X = X.';
end
if (size(Y,2)>size(Y,1))
   Y = Y.';
end

Inn = find(~isnan(X) & ~isnan(Y));
X   = X(Inn);
Y   = Y(Inn);

Corr = corr(X,Y,varargin{:});
N = length(X);

PermInd = ceil(rand(N,Nsim).*N);

switch lower(PermType)
 case 'y'
    SimCorr = corr(X,Y(PermInd),varargin{:});
 case 'x'
    SimCorr = corr(X(PermInd),Y,varargin{:});
 otherwise
    error('Unknown PermInd option');
end


Ngt = length(find(SimCorr>Corr));
Prob = Ngt./Nsim;
