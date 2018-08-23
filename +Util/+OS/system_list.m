function [Stat,Res]=system_list(CommandStr,Files)
% Run the system command on a list of files.
% Package: Util.OS
% Description: Run the system command for a list of files.
% Input  : - System command string (e.g., 'gzip -d %s').
%          - List of files on which to run the system command.
%            This can be any valid input to create_list.m which generate
%            a cell array of files.
% Output : - Vector of Status per file.
%          - Cell array of Results per file.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Stat,Res]=Util.OS.system_list('gzip -d %s',{'File1.gz','File2.gz'});
% Reliable: 2
%--------------------------------------------------------------------------


[~,List] = Util.files.create_list(Files,NaN);
Nlist = numel(List);
Stat  = zeros(Nlist,1);
Res   = cell(Nlist,1);
for Ilist=1:1:Nlist
    [Stat(Ilist),Res{Ilist}] = system(sprintf(CommandStr,List{Ilist}));
end
