function Path=vizquery_path(Prog)
%--------------------------------------------------------------------------
% vizquery_path function                                         Catalogue
% Description: Return the path of the vizquery/cdsclient directory
%              Edit this program before using wget_2mass.m, wget_usnob1.m
%              and wget_ucac4.m.
% Input  : - Optional program name in the cdsclient directory. If provided
%            will check that program exist.
% Output : - Path of the vizquery/cdsclient.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Feb 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Path=vizquery_path;
% Reliable: 1
%--------------------------------------------------------------------------

if (nargin==0),
    Prog = '';
end

% EDIT IF NEEDED
% ProgPath should contain the position of the cdsclient directory
% relative to the position of the vizquery_path.m function.
% Or be empty, if the cdsclient is in the search path
%ProgPath = '../bin/vizquery/cdsclient-3.4/';
ProgPath = '../bin/vizquery/cdsclient-3.80/';

if (isempty(ProgPath)),
    Dir     = '';
    FileSep = '';
else
    Dir = Util.files.which_dir('vizquery_path.m');
    FileSep = filesep;
end
Path     = sprintf('%s%s%s',Dir,FileSep,ProgPath);

if (~isempty(Prog)),
    ProgPath = sprintf('%s%s',Path,Prog);
    if (exist(ProgPath,'file')==0),
        fprintf('Can not find the %s program in the path\n',Prog);
        fprintf('You need to:\n');
        fprintf('  1. install cdsclient: http://cdsarc.u-strasbg.fr/doc/cdsclient.html\n');
        fprintf('  2. Edit vizquery_path.m and edit the ProgPath variable\n');
        error('See vizquery_path.m for instructions');
    end
end
