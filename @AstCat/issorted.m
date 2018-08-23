function [Ans,AstC]=issorted(AstC,Col)
% Check if an AstCat object is sorted by a given column.
% Package: @AstCat
% Description: Check if an AstCat object is sorted by a given column.
% Input  : - AstCat object.
%          - Column name or index to check.
% Output : - Vector of logicals per each AstCat element.
%            True if sorted by column, false otherwise.
%          - AstCat object in which the SortedBy and
%            SortedByCol fields are updated.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Ans,AstC]=issorted(AstC,'YWIN_IMAGE');
% Reliable: 2
%--------------------------------------------------------------------------

            
if (ischar(Col)),
    Col = colname2ind(AstC,Col);
end

Nc  = numel(AstC);
Ans = zeros(size(AstC)); 
for Ic=1:1:Nc,
    Ans(Ic) = issorted(AstC(Ic).Cat(:,Col));
    if (nargout>1),
        if (Ans(Ic)),
            % sorted
            AstC(Ic).SortedBy    = ind2colname(AstC(Ic),Col);
            AstC(Ic).SortedByCol = colname2ind(AstC(Ic),Col);
        else
            AstC(Ic).SortedBy    = [];
            AstC(Ic).SortedByCol = [];
        end
    end
end
