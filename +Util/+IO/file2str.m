function Output=file2str(File,ReadType)
% Read the content of a file into a string or cell vector.
% Package: Util.IO
% Description: Read the content of a file into a string or cell vector
%              (line per element).
% Input  : - File name, cell vector of file names, or a vector of
%            file identifiers.
%            If file identifiers are given, the files will not be closed
%            at the end.
%            If cell vector is given, then each file in each cell will
%            be concatenated.
%          - Read type:
%            'str'  - Read the file including the line breaks, into a
%                     single string (default).
%            'cell' - Read each line in the file (without the line
%                     break), into a cell vector (line per element).
% Output : - String or cell vector containing file content.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 1
%------------------------------------------------------------------------------
DefReadType   = 'str';
if (nargin==1)
   ReadType   = DefReadType;
elseif (nargin==2)
   % do nothing
else
   error('Illegal number of input arguments');
end

% convert 'File' to cell vector:
if (ischar(File)==1)
   CellFile{1} = File;
end
if (iscell(File)==1)
   CellFile    = File;
end

%------------------
%--- open files ---
%------------------
if (isnumeric(File)==1)
   % FID is given
   FID = File;
elseif (iscell(CellFile)==1)
   % get FID for each File
   FID = zeros(length(CellFile),1);
   for If=1:1:length(CellFile)
      FID(If) = fopen(CellFile{If},'r');
   end
else
   error('Unknwon format for File argument');
end

Nf = length(FID);

%------------------
%--- Read files ---
%------------------
switch lower(ReadType)
 case 'str'
    Output = '';
 otherwise
    % do nothing
end

% for each file
K = 0;    % cell line index
for If=1:1:Nf
   switch lower(ReadType)
    case 'str'
       while (feof(FID(If))==0)
          Line = fgets(FID(If));
          Output = sprintf('%s%s',Output,Line);
       end
    case 'cell'
       while (feof(FID(If))==0)
          Line = fgetl(FID(If));
          K = K + 1;
          Output{K} = Line;
       end
    otherwise
       error('Unkown ReadType option');
   end
end


%-------------------
%--- close files ---
%-------------------
if (ischar(File))
   for If=1:1:Nf
      fclose(FID(If));
   end
end

