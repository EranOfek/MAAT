function [Mat,ColCell]=read_csv_with_header(File,varargin)
% Read SDSS CasJobs table into a matrix or table.
% Package: VO.Util
% Description: Read SDSS CasJobs table into a matrix or table.
%              If output is a matrix and some columns contains strings than
%              the unique strings will be replaced by numbers and a
%              dictionary to translate numbers to strings is provided.
% Input  : - File name containing the SDSS CasJobs table.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'StartKey' - A keyword that indicate the end of the header.
%                         If empty (default), then no header.
%            'OutType' - Options are:
%                        'mat' - matrix output. Default.
%                        'table' - table output.
%            'EmptyVal'- Empty val to replace with NaNs. Default is 'null'.
%            'RemDuplicate' - Remove duplicate columns. Default is true.
% Output : - Matrix of table of data.
%          - Cell array of column names.
%          - Column dictionary.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Mat,ColCell]=VO.Util.read_csv_with_header('g13_maged_19me15_1m_20130301_20130331.csv','OutType','table','Delimiter',',','StartKey','data:');
% Reliable: 
%--------------------------------------------------------------------------


DefV.StartKey             = [];
DefV.EndOfLine            = false;
DefV.OutType              = 'mat'; % 'mat' | 'table'
DefV.EmptyVal             = 'null';
DefV.RemDuplicate         = true;
DefV.Delimiter            = ',';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


FID = fopen(File,'r');
Counter = 0;
Line = fgetl(FID);
while isempty(strfind(Line,InPar.StartKey))
    Line = fgetl(FID);
    Counter = Counter + 1;
end
Line = fgetl(FID);
Counter = Counter + 1;


% assume current line contains column names
ColCell = regexp(Line,sprintf('\\%s',InPar.Delimiter),'split');
Ncol    = numel(ColCell);
Format  = '';
for Icol=1:1:Ncol
    Format = sprintf('%s %%s',Format);
end
if (InPar.EndOfLine)
    Format = sprintf('%s\n',Format);
else
    Format = sprintf('%s',Format);
end

NewFile = sprintf('%s_1',File);
system(sprintf('sed "s/%s/%s/g" %s > %s',InPar.EmptyVal,'NaN',File,NewFile));

FID = fopen(NewFile,'r');
C   = textscan(FID,Format,'HeaderLines',Counter+1,'Delimiter',InPar.Delimiter); %,'EmptyValue','null');
fclose(FID);


for Icol=1:1:Ncol
    %C{Icol} = regexprep(C{Icol},InPar.EmptyVal,'NaN');
    Val = str2double(C{Icol});
    if (all(isnan(Val)))
        % do nothing
    else
        C{Icol} = str2double(C{Icol});
    end
end

Mat = table(C{:});

% remove duplicate columns
if (InPar.RemDuplicate)
    [ColCell,UniqueCol] = unique(ColCell);
    Mat = Mat(:,UniqueCol);
end
Mat.Properties.VariableNames = ColCell;
       