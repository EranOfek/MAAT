function [AstC,Col]=colcell2col(AstC)
% Repopulate the AstCat object Col field based on the ColCell field
% Package: @AstCat
% Description: Repopulate the Col field based on the ColCell field.
% Input  : - AstCat object.
% Output : - AstCat object with the Col field repopulated.
%          - The Col field of the last element in the
%            AstCat array.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [AstC,Col]=colcell2col(AstC)
% Reliable: 2
%--------------------------------------------------------------------------

            
Nc = numel(AstC);
for Ic=1:1:Nc
    if (~isempty(AstC(Ic).ColCell))
        AstC(Ic).Col = cell2struct(num2cell(1:1:numel(AstC(Ic).ColCell)),AstC(Ic).ColCell,2);
    end
    Col = AstC(Ic).Col;
end

