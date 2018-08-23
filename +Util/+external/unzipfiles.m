function unzipfiles(zipFilename,files ,outputDirectory)
%   UNZIPFILES Extract the contents of a zip file.
%
%   UNZIPFILES(ZIPFILE) 
%       extracts the contents of a zip file into the current 
%      directory.
%
%   UNZIPFILES(ZIPFILE,FILES) 
%       extracts the files specified in FILES, a string or a
%       a cell array of strings (one or more), from the zipped file
%       into the current directory.
%
%   UNZIPFILES(ZIPFILE,FILES,OUTPUTDIR)  
%       extracts the files specified in FILES, a string or a
%       a cell array of strings, from the zipped file
%       into the directory OUTPUTDIR
%
%   Examples:
%           unzipfiles('junk.zip')
%               will extract the entire contents of junk into the current 
%              current directory.
%          unzipfiles('junk.zip','temp1.mat') 
%               will extract "temp1.mat" into the current 
%              current directory.
%          unzipfiles('junk.zip','temp1.mat','c:\temp') 
%               will extract "temp1.mat" into c:\temp 
%          unzipfiles('junk.zip',[],'c:\temp') 
%               will extract the entire contents of junk into c:\temp
%
%   See also ZIP, UNZIP.


%   This file is an extention of  The MathWorks, Inc's
%   UNZIP function
%    
%   Author :    Tal Pasi
%   Date:     15-Oct-2004


import java.io.*;
import java.util.zip.ZipFile;
import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;

% Argument parsing.
error(nargchk(1,3,nargin));
if (nargin == 1) || (nargin == 2) 
    outputDirectory = pwd;
elseif (nargin == 3) && ~exist(outputDirectory,'dir')
    error('Directory "%s" does not exist.',outputDirectory)
end

if (nargin == 1)
    allFiles =1;
else
    allFiles = 0;
    if isempty(files)
        allFiles = 1;
    else
        % number of files to unzip
        nFiles = size(files,1);
    end
end





% Open the Zip file.
if ~exist(zipFilename,'file')
    error('File "%s" does not exist.',zipFilename);
end
try
    zipFile = ZipFile(zipFilename);
catch
    error('Error opening zip file "%s".',zipFilename);
end

% This InterruptibleStreamCopier is unsupported and may change without notice.
interruptibleStreamCopier = ...
    InterruptibleStreamCopier.getInterruptibleStreamCopier;

% Inflate all entries.
enumeration = zipFile.entries;
while enumeration.hasMoreElements
    zipEntry = enumeration.nextElement;
    % Open output stream.
    if allFiles
        outputName = fullfile(outputDirectory,char(zipEntry.getName));
    else
        % find the wanted files
        for index = 1: nFiles
            if strcmp(zipEntry.getName , files{index})
                outputName = fullfile(outputDirectory,char(files{index}));
                break;
            end
        end
    end
    file = java.io.File(outputName);
    parent = File(file.getParent);
    parent.mkdirs;
    try
        fileOutputStream = java.io.FileOutputStream(file);
    catch
        error('Could not create "%s".',outputName);
    end
    % Extract entry to output stream.
    inputStream = zipFile.getInputStream(zipEntry);
    interruptibleStreamCopier.copyStream(inputStream,fileOutputStream);
    % Close streams.
    fileOutputStream.close;
    inputStream.close;
    
    end
    
end
% Close zip.
zipFile.close;
