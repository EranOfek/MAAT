function S=std_w(X,Sig,Dim)
%--------------------------------------------------------------------------
% std_w function                                                 AstroStat
% Description: Weighted standard deviation.
% Input  : - Matrix of values for which to calculate the weighted std.
%          - Matrix of sigmas (Weight= 1./Sig^2).
%          - Dimension along to calculate the weighted std. Default is 1.
% Output : - Weighted std.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jan 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

if (nargin==2),
    Dim = 1;
end

W = 1./(Sig.^2);  % weight
N = size(X,Dim);

S = sqrt(sum(W.*(X - mean(X,Dim)).^2,Dim)./(((N-1)./N).*sum(W,Dim)));
