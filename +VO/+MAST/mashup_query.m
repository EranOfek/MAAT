function []=mashup_query(varargin)
% SHORT DESCRIPTION HERE
% Package: VO.MAST
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

MashupURL = 'https://mast.stsci.edu/api/v0/invoke?request='; %(Mashup Request Object). 

if nargin==0
    Query = '{''service'':''MAST.Caom.Cone'',''params'':{''ra'':254.28746,''dec'':-4.09933,''radius'':0.2},''format'':''json'',''pagesize'':2000,''removenullcolumns'':True,''timeout'':30,''removecache'':True}';
end
'{'service':'Mast.Caom.Cone','params':{'ra':254.28746,'dec':-4.09933,'radius':0.2},'format':'json','pagesize':2000,'removenullcolumns':True,'timeout':30, 'removecache':True}'

%DefV. = 
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

[T,Status] = urlread(sprintf('%s%s',MashupURL,urlencode(Query)));
