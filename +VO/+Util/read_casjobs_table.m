function [Mat,ColCell,ColDic]=read_casjobs_table(File,varargin)
% Read SDSS CasJobs table into a matrix or table.
% Package: VO.Util
% Description: Read SDSS CasJobs table into a matrix or table.
%              If output is a matrix and some columns contains strings than
%              the unique strings will be replaced by numbers and a
%              dictionary to translate numbers to strings is provided.
% Input  : - File name containing the SDSS CasJobs table.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'OutType' - Options are:
%                        'mat' - matrix output. Default.
%                        'table' - table output.
%            'EmptyVal'- Empty val to replace with NaNs. Default is 'null'.
%            'RemDuplicate' - Remove duplicate columns. Default is true.
% Output : - Matrix of table of data.
%          - Cell array of column names.
%          - Column dictionary.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Mat=VO.Util.read_casjobs_table('sdss_dr14_spec.csv','OutType','table');
% Reliable: 
%--------------------------------------------------------------------------



DefV.OutType              = 'mat'; % 'mat' | 'table'
DefV.EmptyVal             = 'null';
DefV.RemDuplicate         = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


FID = fopen(File,'r');
Line = fgetl(FID);
fclose(FID);

ColCell = regexp(Line,',','split');
Ncol    = numel(ColCell);
Format  = '';
for Icol=1:1:Ncol
    Format = sprintf('%s %%s',Format);
end
Format = sprintf('%s\n',Format);

NewFile = sprintf('%s_1',File);
system(sprintf('sed "s/%s/%s/g" %s > %s',InPar.EmptyVal,'NaN',File,NewFile));

FID = fopen(NewFile,'r');
C   = textscan(FID,Format,'HeaderLines',1,'Delimiter',','); %,'EmptyValue','null');
fclose(FID);

ColDic = Util.struct.struct_def({'Col','UniqueStr'},1,0);

switch lower(InPar.OutType)
    case 'mat'
        Mat = nan(numel(C{1}),Ncol);
        ColI = 0;
        for Icol=1:1:Ncol
            %C{Icol} = regexprep(C{Icol},InPar.EmptyVal,'NaN');
            Val = str2double(C{Icol});
            if (all(isnan(Val)))
                % a string
                ColI = ColI + 1;
                ColDic(ColI).Col       = ColI;
                ColDic(ColI).UniqueStr = unique(C{Icol});
                Nus = numel(ColDic(ColI).UniqueStr);
                for Ius=1:1:Nus
                    C{Icol} = regexprep(C{Icol},ColDic(ColI).UniqueStr{Ius},sprintf('%d',Ius));
                end
            end
            Mat(:,Icol) = str2double(C{Icol});
        end
        % remove duplicate columns
        if (InPar.RemDuplicate)
            [ColCell,UniqueCol] = unique(ColCell);
            Mat = Mat(:,UniqueCol);
        end
    case 'table'
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
    otherwise
        error('Unknown OutType option');
end
    
    
   