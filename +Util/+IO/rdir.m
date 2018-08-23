function Out=rdir(FileTemp)
% recursive dir function
% Package: Util.IO
% Description: run the dir function recursively on a tree of directories
%              looking for specific files.
% Input  : - File template (e.g., '*.mat').
% Output : - Structure array of files information.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Util.IO.rdir('*.mat');
% Reliable: 2
%--------------------------------------------------------------------------


%Out = Util.struct.struct_def({'name','folder','date','bytes','isdir','datenum'},0,1);
D=dir;
Idir=find([D.isdir] & ~strcmp({D.name},'.') & ~strcmp({D.name},'..'));
D=D(Idir);

Nd= numel(D);
Out = dir(FileTemp);
for Id=1:1:Nd
    
    if (D(Id).isdir)
        cd(D(Id).name);
        Out = [Out; Util.IO.rdir(FileTemp)];
        cd ..
    end
end

