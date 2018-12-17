function install_cats(varargin)
% Install the data/+cats catalog directory
% Package: VO.prep
% Description: Install the data/+cats catalog directory. Download
% Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Location' - 
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Dec 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------



DefV.Location             = '~/matlab/data/';
DefV.Dir                  = '+cats';
DefV.URL                  = '';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

PWD = pwd;

cd(InPar.Location);

if exist(InPar.Dir)==7
    % folder exist
else
    mkdir(InPar.Dir);
end

www.pwget(InPar.URL);
T = regexp(InPar.URL,filesep,'split');
untar(T{end});

cd(PWD);
VO.prep.prep_data_dir