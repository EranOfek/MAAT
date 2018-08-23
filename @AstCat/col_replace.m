function AstC=col_replace(AstC,Vec,ColInd,ColName)
% Recplace a column in an AstCat object.
% Package: @AstCat
% Description: Recplace a column in an AstCat object.
% Input  : - AstCat object.
%          - Vector or a matrix or a table columns
%            to replace.
%          - Column indices of the columns to replace.
%          - Cell array of new column names, or a column name.
%            If empty then use the old
%            column name. Default is empty.
% Output : - AstCat object with the columns replaced.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: A=col_replace(A,ones(10,1),2,'new')
% Reliable: 2
%--------------------------------------------------------------------------


           
if (nargin<4),
    ColName = [];
end

if (ischar(ColName)),
    ColName = {ColName};
end

Nc = numel(AstC);
for Ic=1:1:Nc,
    % for each catalog

    ColInd = colname2ind(AstC(Ic),ColInd);
    if (istable(AstC(Ic).Cat)),
        % table treatment
        AstC(Ic).Cat(:,ColInd) = array2table(Vec);
    else
        AstC(Ic).Cat(:,ColInd) = Vec;
    end

    if (~isempty(ColName) && ~isempty(AstC(Ic).ColCell)), 
        % replace column names
        AstC(Ic).ColCell(ColInd) = ColName(:);
        AstC(Ic) = colcell2col(AstC);
    end

end
            
