function [Corr,Prob,Num]=corrsim_cov(Mat,Nsim,PermType,varargin)
% Correlation matrix between N columns and bootstrap estimation of errors.
% Package: Util.stat
% Description: Given a matrix with N columns, calculate the correlation
%              between each pair of columns and use the
%              bootstrap method to estimate the probability to get
%              a correlation larger than the observed correlation.
%              The function ignores NaN values.
% Input  : - Matrix with N columns.
%          - Number of bootstrap simulations. Default is 1000.
%          - Which vector to permute in bootstrap test {'x'|'y'}.
%            Default is 'y'.
%          * Arbitrary number of input argument to pass to the
%            built in function corr.m
%            'type' - 'Pearson' (the default) to compute Pearson's linear
%                     correlation coefficient.
%                     'Kendall' to compute Kendall's tau.
%                     'Spearman' to compute Spearman's rho.
% Output : - Matrix of all correlations.
%          - Matrix of all probabilities to get larger correlation.
%          - Matrix of number of meauseremnets used in each correlation.
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Mar 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Corr,Prob]=corrsim_cov(rand(100,4));
% Reliable: 1
%------------------------------------------------------------------------------

Def.Nsim     = 1000;
Def.PermType = 'Y';

if (nargin==1)
   Nsim  = Def.Nsim;
   PermType = Def.PermType;
elseif (nargin==2)
   PermType = Def.PermType;
else
   % do nothing
end


[M,N] = size(Mat);
Corr  = zeros(N,N).*NaN;
Prob  = zeros(N,N).*NaN;

for I=1:1:N
    for J=1:1:I-1
        [C,P,SimCorr,N]=corrsim(Mat(:,I),Mat(:,J),Nsim,PermType,varargin{:});
        Corr(I,J) = C;
        Prob(I,J) = P;
        Num(I,J)  = N;
    end
end


