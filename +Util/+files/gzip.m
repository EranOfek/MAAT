function OutFile=gzip(InFile,Args)
% Execute gzip on files, or a cell array of files
% Package: AstroUtil.files
% Description: Execute the gzip system command (file compression) on files,
%              or a list of files.
% Input  : - String containing a file name, file name with wild cards or
%            ranges, a file (starting with @) containing a list of files,
%            or a cell array of file names. See create_list for options.
%          - Additional arguments to pass to the gzip command.
%            Default is ''.
% Output : - Cell array of compressed file names (i.e., '.gz' extension
%            added).
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
    system(sprintf('gzip %s %s',Args,OutFile{Ifile}));
    OutFile{Ifile} = sprintf('%s.gz',OutFile{Ifile});
end


