function fit_roblin(X,Y)
%--------------------------------------------------------------------------
% fit_roblin function                                               FitFun
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


% NumVarargs = length(varargin);
% if NumVarargs > 3
%      errId = 'fit_roblin.m:TooManyInputArguments';
%      errMsg = 'InPar1, [InPar2, InPar3]';
%      error(errId, errMsg);
% end
% Gaps = cellfun(@isempty, varargin);
% DefArgs = {InPar1Def InPar2Def InPar3Def};    % default input arguments
% Suboptargs = DefArgs(1 : NumVarargs);
% varargin(Gaps) = Suboptargs(Gaps);
% DefArgs(1 : NumVarargs) = varargin;
% [Par1 Par2 Par3] = DefArgs{:}
% 
% 
% DefV. = 
% InPar = set_varargin_keyval(DefV,'n','use',varargin{:});


Flag = X>0 & Y>0;
LogX = log10(X(Flag));
LogY = log10(Y(Flag));

[MaxLX] = max(LogX);

B = binning([LogX,LogY],0.5,-1,ceil(MaxLX));
D = diff(10.^(B(:,[4 8])))
plot(10.^B(:,4),10.^B(:,8),'ro')

log10(D(:,2)./D(:,1))

