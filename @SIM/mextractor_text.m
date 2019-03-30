function []=mextractor_text(varargin)
% SHORT DESCRIPTION HERE
% Package: @SIM
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


NumVarargs = length(varargin);
if NumVarargs > 3
     errId = 'mextractor_text.m:TooManyInputArguments';
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
InPar = InArg.populate_keyval(DefV,varargin,mfilename);
