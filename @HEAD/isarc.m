function []=isarc(varargin)
%--------------------------------------------------------------------------
% isarc function                                               class/@HEAD
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
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


NumVarargs = length(varargin);
if NumVarargs > 3
     errId = 'isarc.m:TooManyInputArguments';
     errMsg = 'InPar1, [InPar2, InPar3]';
     error(errId, errMsg);
end
Gaps = cellfun(@isempty, varargin);
DefArgs = {InPar1Def InPar2Def InPar3Def};    % default input arguments
Suboptargs = DefArgs(1 : NumVarargs);
varargin(Gaps) = Suboptargs(Gaps);
DefArgs(1 : NumVarargs) = varargin;
[Par1 Par2 Par3] = DefArgs{:}


DefV. = 
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
