function Dir=dir_cell(Cell)
% dir like command for a cell of file names.
% Package: Util.files
% Description: dir like command for a cell of file names.
% Input  : - A cell array of file names.
% Output : - Dir output.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Sep 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: A=Util.files.dir_cell({'GI5_048014_GAMA12_13209-nd-rrhr.fits','GI5_048014_GAMA12_13209-nd-cnt.fits'})
% Reliable: 2
%--------------------------------------------------------------------------

Ncell = numel(Cell);
Dir   = Util.struct.struct_def({'name','folder','date','bytes','isdir','datenum'},Ncell,1);
%Dir   = Util.struct.struct_def({'name','date','bytes','isdir','datenum'},Ncell,1);

for Icell=1:1:Ncell
    if (exist(Cell{Icell},'file')>0 || exist(Cell{Icell},'dir')>0)
        dir(Cell{Icell})
        Dir(Icell) = dir(Cell{Icell});        
    end
end