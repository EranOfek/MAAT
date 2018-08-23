function copy_files_from_dirtree(Location,CopyTo,CopyType)
% Copy or movde all files recursively in a directory tree.
% Package: Util.files
% Description: Given a location (or the present working directory), look
%              for all files in the subdirectories and copy them to the main
%              directory.
% Input  : - Parent directory. Default is the present working directory.
%          - Directory to copy the data to. Default is the present working
%            directory.
%          - CopyType: 'move' - use move instead of copy.
%                     default is 'copy'.
% Output : null
% Tested : Matlab 7.13
%     By : Eran O. Ofek                    Jun 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 1
%--------------------------------------------------------------------------

Def.CopyType = 'copy';
if (nargin==0),
   Location = pwd;
   CopyTo   = pwd;
   CopyType = Def.CopyType;
elseif (nargin==1),
   CopyTo   = pwd;
   CopyType = Def.CopyType;
elseif (nargin==2),
   CopyType = Def.CopyType;
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end

PWD = pwd;
cd(Location);

Files = dir('*');
Nf    = length(Files);
for If=3:1:Nf,
   switch (Files(If).isdir)
    case 0
       % copy file
       switch lower(CopyType)
        case {'copy','cp'}
           copyfile(Files(If).name,CopyTo);
        case {'move','mv'}
           movefile(Files(If).name,CopyTo);
        otherwise
	   error('Unknown CopyType option');
       end
    case 1
       % call function recursively
       cd(Files(If).name);
       copy_files_from_dirtree(pwd,CopyTo,CopyType);
    otherwise
       error('Unknown isdir option');
   end
end

% go back
%cd ..
cd(PWD);