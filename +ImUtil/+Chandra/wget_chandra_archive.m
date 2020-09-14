function []=wget_chandra_archive(varargin)
% SHORT DESCRIPTION HERE
% Package: ImUtil
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV.BaseURL              = 'ftp://cda.cfa.harvard.edu/pub/cal';
DefV.AO                   = 'all';  % or number
DefV.BaseDir              = pwd;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

LenBase = numel(InPar.BaseURL);

URL = sprintf('%s/ao%02d/',InPar.BaseURL,InPar.AO);
%F=www.ftp_dir_list('ftp://cda.cfa.harvard.edu/pub/cal/',false)
%F = www.ftp_dir_list('ftp://cda.cfa.harvard.edu/pub/cal/ao07/cat9/',true);
F = www.ftp_dir_list(URL,true);

% move to directory in which to save the data
PWD = pwd;
cd(InPar.BaseDir);


Nf = numel(F);
for If=1:1:Nf
    
    FileName       = F(If).name;
    FileNameAndDir = F(If).URL(LenBase+1:end);
    Tmp = regexp(F(If).URL,'/','split');
    ObsID = str2double(Tmp{8});
    CatN  = Tmp{7};
    AO    = Tmp{6};
    
    www.pwget(F(If).URL)
end

% go back to home dir
cd(PWD);
