function File=load_from_zip(ZipFile,List)
% Extract and load files from a zip file.
% Package: Util.IO
% Description: Extract and load files from a zip file.
% Input  : - Zip file name.
%          - A cell array of files to extract and load from the zip file.
% Output : - A cellay array of the loaded files.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    May 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

String = '';
for I=1:1:length(List)
    String = sprintf('%s %s',String,List{I});
end
TmpDir = '/tmp';

system(sprintf('unzip %s %s -d %s',ZipFile,String,TmpDir))

File = cell(numel(List),1);
for I=1:1:length(List)
    File{I} = load2(sprintf('%s%s%s',TmpDir,filesep,List{I}));
    delete(sprintf('%s%s%s',TmpDir,filesep,List{I}));
end
