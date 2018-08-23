function OutFile=gunzip(InFile,Args)
% Execute gunzip on files, or a cell array of files
% Package: AstroUtil.files
% Description: Execute the gunzip system command (file uncompression) on
%              files, or a list of files.
% Input  : - String containing a file name, file name with wild cards or
%            ranges, a file (starting with @) containing a list of files,
%            or a cell array of file names. See create_list for options.
%          - Additional arguments to pass to the gunzip command.
%            Default is ''.
% Output : - Cell array of compressed file names (i.e., '.gz' extension
%            removed).
%     By : Eran O. Ofek                    Oct 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstroUtil.files.gzip('*.fits');
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1),
    Args = '';
end

[~,OutFile] = Util.files.create_list(InFile,NaN);
Nfile = numel(OutFile);
for Ifile=1:1:Nfile,
    system(sprintf('gunzip %s %s',Args,OutFile{Ifile}));
    OutFile{Ifile} = regexprep(OutFile{Ifile},'.gz','');
end


