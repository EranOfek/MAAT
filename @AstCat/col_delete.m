function AstC=col_delete(AstC,Col)
% Delete specific columns from an AstCat object.
% Package: @AstCat
% Description: Delete specific columns from an AstCat class object.
% Input  : - AstCat class object.
%          - Columns indices or a cell array of column
%            names to delete.
% Output : - AstCat class object with the deleted columns.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstC=col_select(AstC,{'XWIN_IMAGE','YWIN_IMAGE'})
% Reliable: 2
%--------------------------------------------------------------------------


CatField         = 'Cat';
%ColField         = 'Col';
ColCellField     = 'ColCell';
ColUnitsField    = 'ColUnits';
SortedByField    = 'SortedBy';
SortedByColField = 'SortedByCol';


Nc = numel(AstC);
for Ic=1:1:Nc,
    ColIndDel = colname2ind(AstC(Ic),Col);
    Ncol      = numel(AstC(Ic).(ColCellField));
    % column indices to keep
    ColInd    = setdiff((1:1:Ncol),ColIndDel);

    AstC(Ic).(CatField) = AstC(Ic).(CatField)(:,ColInd);
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
