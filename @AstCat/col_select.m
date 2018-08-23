function AstC=col_select(AstC,Col,Flag)
% Select specific columns from an AstCat class object.
% Package: @AstCat
% Description: Select specific columns from an AstCat class object.
% Input  : - AstCat class object.
%          - Columns indices or a cell array of column
%            names to select. Empty will return no columns.
%            If Inf then will return all columns.
%          - Optional logical flags or indices of lines to return.
%            If empty, then return all lines.
%            If a number in a string then will return the last N rows.
%            Default is empty.
% Output : - AstCat class object with the selected columns.
% See also: AstCat/show.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstC=col_select(AstC,Col)
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<3),
    Flag = [];
end

CatField         = 'Cat';
%ColField         = 'Col';
ColCellField     = 'ColCell';
ColUnitsField    = 'ColUnits';
SortedByField    = 'SortedBy';
SortedByColField = 'SortedByCol';

Nc = numel(AstC);
for Ic=1:1:Nc,
    
    if (isempty(Flag)),
        FlagIc = true(size(AstC(Ic).(CatField),1),1);
    
    else
        if (ischar(Flag)),
            Nline   = size(AstC(Ic).(CatField),1);
            FlagIc  = str2double(Flag);
            FlagIc  = (max(Nline-FlagIc,1):Nline)';
        else
            FlagIc = Flag;
        end
    end
    
    if (iscell(Col)),
        ColInd = colname2ind(AstC(Ic),Col);
    else
        if (isinf(Col)),
            % Col is infinity - return all columns
            ColInd = (1:1:size(AstC(Ic).(CatField),2));
        else
            ColInd = colname2ind(AstC(Ic),Col);
        end
    end
    
    AstC(Ic).(CatField) = AstC(Ic).(CatField)(FlagIc,ColInd);
    if (~isempty(AstC(Ic).(ColCellField))),
        AstC(Ic).(ColCellField) = AstC(Ic).(ColCellField)(ColInd);
        AstC(Ic) = colcell2col(AstC(Ic));
    end
    if (~isempty(AstC(Ic).(ColUnitsField))),
        AstC(Ic).(ColUnitsField) = AstC(Ic).(ColUnitsField)(ColInd);
    end
    AstC(Ic).(SortedByField) = [];
    AstC(Ic).(SortedByColField) = [];
end
