function [AstC]=astcat_array2table(AstC)
%--------------------------------------------------------------------------
% astcat_array2table function                                class/@AstCat
% Description: Convert the catalog field in an AstCat object to a table.
% Input  : - AstCat object.
% Output : - AstCat object in which the Cat field is a table.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstC=astcat_array2table(AstC);
% Reliable: 2
%--------------------------------------------------------------------------

Nc = numel(AstC);
for Ic=1:1:Nc,
    if (isempty(AstC(Ic).ColCell)),
        AstC(Ic).Cat = array2table(AstC(Ic).Cat);
    else
        AstC(Ic).Cat = array2table(AstC(Ic).Cat,'VariableNames',AstC(Ic).ColCell);
    end
end
    