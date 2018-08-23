function [URLs]=url2url_key_val(URLpar,varargin)
% SHORT DESCRIPTION HERE
% Package: www
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV.Amp                  = '&amp;';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            
if (~iscell(URLpar))
    URLpar = {URLpar};
end

Nurl = numel(URLpar);
for Iurl=1:1:Nurl
    SplitURL = regexp(URLpar{Iurl},'?','split');
    URLs(Iurl).URL = SplitURL{1};
    SplitPar = regexp(SplitURL{2},InPar.Amp,'split');
    Npar = numel(SplitPar);
    KeyVal = cell(0,1);
    for Ipar=1:1:Npar
        SplitKV = regexp(SplitPar{Ipar},'=','split');
        KeyVal  = [KeyVal, SplitKV{:}];
    end
        
    URLs(Iurl).KeyVal = KeyVal;
end

    
