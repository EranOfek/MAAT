function AstC=col_duplicate(AstC,Cols,NewColNames)
%--------------------------------------------------------------------------
% col_duplicate function                                     class/@AstCat
% Description: Duplicate specific columns in an AstCat object into the
%              additional columns at the end of the catalog.
% Input  : - An AstCat object.
%          - Cell array of column names or vector of column indices to
%            duplicate.
%          - Cell array of names for the duplicated columns.
%            Defualt is "dup_(OriginalName)".
% Output : - An AstCat object with the duplicated columns.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstC=col_duplicate(AstC,{'XWIN_IMAGE','YWIN_IMAGE'});
% Reliable: 
%--------------------------------------------------------------------------

Prefix = 'dup_';
if (nargin<3)
    NewColNames = [];
end

CatField     = 'Cat';

Nc = numel(AstC);
for Ic=1:1:Nc
    % for each catalog
    ColInd = colname2ind(AstC(Ic),Cols);
    
    if (isempty(NewColNames))
        NewColNamesCell = Util.cell.cellstr_prefix(AstC(Ic).ColCell(ColInd),Prefix);
    else
        NewColNamesCell = NewColNames;
    end
    
    AstC(Ic).(CatField) = [AstC(Ic).(CatField), AstC(Ic).(CatField)(:,ColInd)];
    AstC(Ic).ColCell    = [AstC(Ic).ColCell(:)', NewColNamesCell];
    AstC                = colcell2col(AstC);
end

    
