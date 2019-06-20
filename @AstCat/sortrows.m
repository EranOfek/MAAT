function [AstC,SI]=sortrows(AstC,Col)
% Sort an AstCat object by a given column.
% Package: @AstCat
% Description: Sort an AstCat object by a given column.
% Input  : - AstCat object.
%          - Column name or index to sort by.
%            If negative column index than sort by descending
%            order.
%            Also can be a cell array of names or a vector.
% Output : - Sorted AstCat object.
%            Also updated the SortedBy and SortedByCol fields.
%          - Indices of the original indices sorted.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstC=sortrows(AstC,[2 -3]);
%          AstC=sortrows(AstC,'YWIN_IMAGE');
%          AstC=sortrows(AstC,{'YWIN_IMAGE','XWIN_IMAGE'});
% Reliable: 2
%--------------------------------------------------------------------------


            

Nc = numel(AstC);
for Ic=1:1:Nc
    if (ischar(Col))
        Col = colname2ind(AstC(Ic),Col);
    end
    if istable(AstC(Ic).Cat)
        AstC(Ic).Cat         = sortrows(AstC(Ic).Cat,Col);
    else
        [AstC(Ic).Cat,SI]    = sortrows(real(AstC(Ic).Cat),Col);
    end
    AstC(Ic).SortedBy    = ind2colname(AstC(Ic),Col);
    AstC(Ic).SortedByCol = colname2ind(AstC(Ic),Col);
end
