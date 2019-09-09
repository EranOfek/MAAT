function [Table]=read_votable(FileStr,varargin)
% XML VO table reader 
% Package: VO.Util
% Description: Simple XML (VO table) reader.
% Input  : - Either a cell containing a string containing of file name,
%            a file identifier, or a string containing the VO table string.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Table]=VO.Util.read_votable({'filename'});
%          [Table]=VO.Util.read_votable(FID);
%          [Table]=VO.Util.read_votable(String);
% Reliable: 2
%--------------------------------------------------------------------------

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if iscellstr(FileStr)
    % file name
    FID = fopen(FileStr{1},'r');
    Str = fscanf(FID,'%s',Inf);
    fclose(FID);
elseif isnumeric(FileStr)
    % FID
    Str = fscanf(FileStr,'%s',Inf);
elseif ischar(FileStr)
    % string
    Str = FileStr;
else
    error('Unknown FileStr format');
end


ColNames = regexp(Str,'<FIELD name="(?<name>\w+)" datatype="(?<datatype>\w+)"','names');
Ncol = numel(ColNames);

Val  = regexp(Str,'<TD>(?<val>[\w\-\+\.]+)</TD>','names');
Nval = numel(Val);
Nline = Nval./Ncol;

if (Nval==0)
    Table = [];
else
    if (Nline~=floor(Nline))
        error('likely problem in reading (e.g., special characters)');
    end
    V = reshape({Val.val},Ncol,Nline).';
    Table = cell2table(V);
    Table.Properties.VariableNames = {ColNames.name};
    for Icol=1:1:Ncol
        Col = Table.Properties.VariableNames{Icol};

        switch lower(ColNames(Icol).datatype)
            case {'double','float'}
                Table.(Col) = str2double(Table.(Col));

            case {'char'}
                % do nothing - already a string
            case {'int','long'}
                Table.(Col) = str2double(Table.(Col));
            case 'unsignedbyte'
                for Iline=1:1:Nline
                    Table.(Col){Iline} = Table.(Col){Iline}(3:end);
                end
                Table.(Col) = hex2dec(Table.(Col));

            otherwise
                fprintf('Need to define a new datatype: %s\n',ColNames(Icol).datatype);
                error('Unknown datatype');
        end
    end

end
