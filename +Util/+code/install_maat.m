function []=install_maat(DirToInstallMMAT)
% Install
% Package: Util.code
% Description:
% Input  : -
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : -
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2017
%     Edited: Ygal Y. Klein 29.5.2017
%    URL : http://weizmann.ac.il/home/eofek/matlabn/
% Example:
% Reliable:
%--------------------------------------------------------------------------

FunLink  = 'https://webhome.weizmann.ac.il/home/eofek/matlab/Fun.tar.gz';
CatsLink = 'https://webhome.weizmann.ac.il/home/eofek/matlab/data_cats.tar.gz';
StartLink= 'https://webhome.weizmann.ac.il/home/eofek/matlab/startup.m';

PWD = pwd;

if (nargin==0)
  DirToInstallMMAT = [];
end

if (isempty(DirToInstallMMAT))
   DirToInstallMMAT = input('Type full-path directory name (whitout commas at head and tail) in which you want to instal MAAT (default is ~/matlab): ','s');
end

if (isempty(DirToInstallMMAT))
    DirToInstallMMAT = '~/matlab';
end

if (~exist(DirToInstallMMAT,'dir'))
    error('Directory %s does not exist - create and start again');
end

cd(DirToInstallMMAT)

fprintf('Downloading the functions tar file\n');
untar(FunLink);
fprintf('Downloading the Data/+cats tar file\n');
untar(CatsLink);

Input = input('Do you want to install the startup.m file in the matlab/ directory Y/N (default is N): ','s');
switch lower(Input)
    case 'y'
        
        [~,name,ext] = fileparts(StartLink);
        StartupFullPath = websave([name ext], StartLink);
        
        
    otherwise
        % do nothing
end

fprintf('MAAT is installed\n');

end
