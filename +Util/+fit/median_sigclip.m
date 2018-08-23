function RM=median_sigclip(X,Perc)
%--------------------------------------------------------------------------
% median_sigclip function                                          General
% Description: A robust median calculation, by removing lower and upper
%              percentiles prior to the median calculation.
% Input  : - Vector of values.
%          - Lower and upper percentiles [low, upper].
%            Default is [0.2 0.8].
% Output : - Robust median.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Exampel: RM=median_sigclip(rand(1000,1));
% Reliable: 2
%--------------------------------------------------------------------------

Def.Perc = [0.2 0.8];
if (nargin==1)
    Perc = Def.Perc;
end

RM = nanmedian(X(X>quantile(X,Perc(1)) & X<quantile(X,Perc(2))));







