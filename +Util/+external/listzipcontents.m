function filelist = listzipcontents(zipFilename)
%LISTZIPCONTENTS Lists the contents of zip file.
%
%   LISTZIPCONTENTS(ZIPFILENAME) lists the archived contents of ZIPFILENAME
%
%   ZIPFILENAME is a string specifying the name of the zip file. 
%   ZIPFILENAME can include the directory name; otherwise, the file must be
%   in the current directory or in a directory on the MATLAB path.
%
%   FILENAMES = LISTZIPCONTENTS(ZIPFILENAME) lists the zip file contents
%   and returns the file names into the string cell array FILENAMES.
%
%   Unsupported zip files
%   ---------------------
%   LISTZIPCONTENTS does not support password-protected or encrypted zip 
%   archives.
%
%   Examples
%   --------
%   % List the contents of demos.zip 
%   listzipcontents('demos.zip')
%
%   B.C. Hamans, UMC St Radboud, 2011 (B.C.Hamans@rad.umcn.nl)
%
%   See also FILEATTRIB, GZIP, GUNZIP, TAR, UNTAR, UNZIP, ZIP.

%Copyright (c) 2011, Bob Hamans
%All rights reserved.

%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:

%* Redistributions of source code must retain the above copyright
%notice, this list of conditions and the following disclaimer.
%* Redistributions in binary form must reproduce the above copyright
%notice, this list of conditions and the following disclaimer in
%the documentation and/or other materials provided with the distribution

%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%                       SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%                       INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error(nargchk(1,1,nargin,'struct'));

% Create a Java file of the ZIP filename.
zipJavaFile  = java.io.File(zipFilename);
% Create a Java ZipFile and validate it.
zipFile = org.apache.tools.zip.ZipFile(zipJavaFile);
% Extract the entries from the ZipFile.
entries = zipFile.getEntries;

% Initialize the file list.
filelist={};

% Loop through the entries and add to the file list.
while entries.hasMoreElements
    filelist = cat(1,filelist,char(entries.nextElement));
end
end

