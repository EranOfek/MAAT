function [List,FileNames]=find_urls_ftp(URL,varargin)
% Find files in a FTP link
% Package: www
% Description: Make a list of files in an FTP link.
% Input  : - FTP URL.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - Cell array of links names.
%          - Cell array of file names.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [List,FileNames]=www.find_urls_ftp(URL);
% Reliable: 2
%--------------------------------------------------------------------------

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Str = urlread(URL);
C   = textscan(Str,'%s %d %d %d %d %s %d %s %s\n');
FileNames = C{end};
Nf        = numel(FileNames);
List      = cell(Nf,1);
for If=1:1:Nf
    List{If} = sprintf('%s%s',URL,FileNames{If});
end

