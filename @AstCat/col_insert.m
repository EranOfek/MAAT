function AstC=col_insert(AstC,Vec,ColInd,ColName)
% Insert a column into an AstCat object.
% Package: @AstCat
% Description: Insert a column into an AstCat object.
% Input  : - AstCat object.
%          - Vector or a cell vector or a table column
%            to insert.
%          - Column index in which to insert the new column.
%          - Column name of the new column.
% Output : - AstCat object with new column.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AA=col_insert(A,ones(10,1).*3,4,'CAA')
% Reliable: 2
%--------------------------------------------------------------------------


           
if (numel(ColInd)>1),
    error('col_insert can insert one column at a time');
end
            
Nc = numel(AstC);
for Ic=1:1:Nc,
    
    if istable(AstC(Ic).Cat),
        AstC(Ic).Cat = [AstC(Ic).Cat(:,[1:ColInd-1]), array2table(Vec), AstC(Ic).Cat(:,[ColInd:end])];
    else
        AstC(Ic).Cat = [AstC(Ic).Cat(:,[1:ColInd-1]),             Vec,  AstC(Ic).Cat(:,[ColInd:end])];
    end
    AstC(Ic).ColCell = [AstC(Ic).ColCell(1:ColInd-1), ColName, AstC(Ic).ColCell(ColInd:end)];
    AstC(Ic)         = colcell2col(AstC(Ic));
end
