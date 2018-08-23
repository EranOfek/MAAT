function AstC=row_delete(AstC,Ind)
% Delete rows from a single element AstCat object.
% Package: @AstCat
% Description: Given a single AstCat object, return a new AstCat with
%              specific rows deleted.
% Input  : - A single AstCat object.
%          - Indices of rows to delete.
% Output : - An AstCat object with the deleted rows.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstC = row_delete(AstC,(1:5));
% Reliable: 
%--------------------------------------------------------------------------


if (numel(AstC)>1),
    error('AstCat object must contain a single catalog');
end

Nrow = size(AstC.Cat,1);
Ind  = setdiff((1:1:Nrow).',Ind);
AstC.Cat = AstC.Cat(Ind,:);
