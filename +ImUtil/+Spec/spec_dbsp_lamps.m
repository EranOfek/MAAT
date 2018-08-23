function ArcName=spec_dbsp_lamps(KeyVal)
%--------------------------------------------------------------------------
% spec_dbsp_lamps function                                        ImSpec
% Description: 
% Input  : - 
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% Tested : Matlab R2011b
%     By : Eran O. Ofek                       Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


KeyVal = Util.string.spacedel(KeyVal);

% retrieve lamp names
ArcKeyValue = {'0000000','1000000','0100000','0010000','0001000','0000100','0000010','0000001','0001100','0001110'};
ArcNameList = {NaN      ,'D',      'Fe+Ar',  'Hg',     'Ar',     'Ne',     'He',     'InCand' ,'Ne+Ar'  ,'Ne+Ar'};

ArcName = ArcNameList{strcmpi(KeyVal,ArcKeyValue)};

