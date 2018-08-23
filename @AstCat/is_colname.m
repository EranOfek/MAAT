function IsCol=is_colname(AstC,ColNames,Case)
% Check if a string or a cell array of string is a valid AstCat column name
% Package: @AstCat
% Description: Check if a string or a cell array of string is a valid
%              column name in an AstCat object.
% Input  : - AstCat object.
%          - A string of column name or a cell array of column names.
%          - Case senstive {true|false}. Default is true.
% Output : - Matrix of logical flags indicating if a column name exist
%            in the AstCat object.
%            Each row corresponds to one AstCat element, while each
%            column corresponds to the column name
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: IsCol=is_colname(AstC,{'X','Y'});
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<3),
    Case = true;
end

if (~iscell(ColNames)),
    ColNames = {ColNames};
end
Ncol = numel(ColNames);

if (Case)
    FunStrCmp = @strcmp;
else
    FunStrCmp = @strcmpi;
end


Ncat  = numel(AstC);
IsCol = false(Ncat,Ncol);
for Icat=1:1:Ncat,
    % for each element
    for Icol=1:1:Ncol,
        % for each requested column
        IsCol(Icat,Icol) = any(FunStrCmp(ColNames{Icol},AstC(Icat).ColCell));
    end
end

    