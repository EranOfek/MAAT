function [AstC]=astcat_table2array(AstC)
%--------------------------------------------------------------------------
% astcat_array2table function                                class/@AstCat
% Description: Convert the catalog field in an AstCat object from a table
%              to array.
% Input  : - AstCat object.
% Output : - AstCat object in which the Cat field is an array.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstC=astcat_table2array(AstC);
% Reliable: 2
%--------------------------------------------------------------------------

Nc = numel(AstC);
for Ic=1:1:Nc,
    if (isempty(AstC(Ic).ColCell)),
        AstC(Ic).ColCell = AstC(Ic).Cat.Properties.VariableNames;
        AstC(Ic).Cat     = table2array(AstC(Ic).Cat);
        AstC(Ic)         = colcell2col(AstC(Ic)); 
    else
        AstC(Ic).Cat = table2array(AstC(Ic).Cat);
    end
end
    