function [Nrow,Ncol]=sizecat(AstC)
% Return the size of the catalog of each AstCat object.
% Package: @AstCat
% Description: Return the size of the catalog of each AstCat object.
% Input  : - AstCat object.
% Output : - Matrix of number of rows, one per AstCat element.
%          - Matrix of number of columns, one per AstCat element.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [N,M] = sizecat(AstC)
% Reliable: 2
%--------------------------------------------------------------------------


Nc = numel(AstC);
Nrow = zeros(size(AstC));
Ncol = zeros(size(AstC));
for Ic=1:1:Nc,
    [Nrow(Ic), Ncol(Ic)] = size(AstC(Ic).Cat);
end