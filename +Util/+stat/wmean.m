function [M,E,S]=wmean(Vec,Err,Dim,IgnoreNaN)
%--------------------------------------------------------------------------
% wmean function                                                 AstroStat
% Description: Calculated the weighted mean of a sample.
% Input  : - Either a two column matrix [Val, Err] or a matrix of values,
%            while the errors are in the second argument.
%          - Optional mtarix of errors. If given, then the first input
%            argument is treated as values.
%          - If the first two input arguments are provided than this is the
%            dimension along to calculate the weighted mean.
%            Default is 1.
%          - Ignore nans. Default is true.
% Output : - Weighted mean.
%          - Weighted error on weighted mean.
%          - Weighted standard deviation.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jun 1998
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [M,E]=Util.stat.wmean([1;2;3;4],[1;1;1;50]);
% Reliable: 2
%--------------------------------------------------------------------------
ColVal  = 1;
ColErr  = 2;

if (nargin<4)
    IgnoreNaN = true;
end

if (nargin>1)
    % Err is given in second argument
    % do nothing
    if (nargin==2)
        Dim = 1;
    end
else
    % Vec is a two coumn vector
    Err = Vec(:,ColErr);
    Vec = Vec(:,ColVal);
    Dim = 1;
end

% Ignore NaNs
if (IgnoreNaN)
    Flag = ~isnan(Vec) & ~isnan(Err);
    Vec  = Vec(Flag);
    Err  = Err(Flag);
end

%E = sqrt(1./sum(1./VecValue(I).^2));
E = sqrt(1./sum(1./Err.^2,Dim));
M = sum(Vec./(Err.^2),Dim)./sum(1./(Err.^2),Dim);
W = 1./Err.^2;  % weight
S = sqrt((sum(W.*Vec.^2).*sum(W) - sum(W.*Vec).^2)./(sum(W).^2 - sum(W.^2)));

