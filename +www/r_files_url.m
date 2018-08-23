function [List,IsDir,FileName]=r_files_url(URL,varargin)
% Recursively get links to all files in www directory list.
% Package: www
% Description: Recursively get links to all files in www directory list.
%              Given a URL that contains a directory tree with files,
%              get the links to all file names in the directory tree.
% Input  : - String containing URL.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            see www.find_urls.m for options.
% Output : - A cell array of links found within the URL.
%          - A flag indicating if link is a possible directory.
%            Directories are identified by the '/' sign in the end
%            of the link.
%          - A cell array of file names in the URL name (i.e., string after
%            last "/").
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Nov 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% [LL,IsDir,FN]=www.r_files_url('https://ztfweb.ipac.caltech.edu/ztf/archive/sci/2017/1015/');
% Reliable: 2
%--------------------------------------------------------------------------

[ListU,IsDirU,FileNameU] = www.find_urls(URL,varargin{:});
Ilink1     = Util.cell.isempty_cell(regexp(FileNameU,'\?C=','match'));
Ilink2     = ~Util.cell.isempty_cell(regexp(ListU,'\.','match'));
Ilink      = find(Ilink1 & Ilink2 & IsDirU);
Ifile      = find(Ilink1 & Ilink2 & ~IsDirU);

% prepare list of files
List       = ListU(Ifile);
IsDir      = IsDirU(Ifile);
FileName   = FileNameU(Ifile);


Ndir   = numel(Ilink);
for Idir=1:1:Ndir
    [ListD,IsDirD,FileNameD] = www.r_files_url(ListU{Ilink(Idir)},varargin{:});
    List           = [List; ListD];
    IsDir          = [IsDir; IsDirD];
    FileName       = [FileName; FileNameD];
end
