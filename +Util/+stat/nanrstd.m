function StD=nanrstd(Data,Flag,Dim)
% Robust nanstd.
% Package: Util.stat
% Description: Robust nanstd. Estimating the std (like nanstd.m), based
%              on the 50-percentile of the distribution.
% Input  : - Data.
%          - If 0, normalize by N-1. If 1, normalize by N. Default is 0.
%          - Dimension along to calculate StD. Default is 1.
% Outout : - Robust Std (output and behavior is like nanstd.m).
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: StD=nanrstd(randn(5,5))
% Reliable: 2
%--------------------------------------------------------------------------

Def.Flag = 0;
Def.Dim  = 1;
if (nargin==1)
    Flag  = Def.Flag;
    Dim   = Def.Dim;
elseif (nargin==2)
    Dim   = Def.Dim;
else
    % do nothing
end

Low=prctile(Data,25,Dim);
Upp=prctile(Data,75,Dim);

StD=(Upp-Low)./1.34897950039216;  % the constant is: norminv(0.75,0,1).*2