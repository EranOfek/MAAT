function [AstC]=astcat_flipud(AstC)
%--------------------------------------------------------------------------
% astcat_flipud function                                     class/@AstCat
% Description: flip up down of each catalog in an AstCat object.
% Input  : - An AstCat object.
% Output : - An AstCat object in which the Cat field is fliped up/down,
%            and the SortedBy fields are set to empty.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [AstC]=astcat_flipud(AstC)
% Reliable: 
%--------------------------------------------------------------------------


Nc = numel(AstC);
for Ic=1:1:Nc,
    AstC(Ic).Cat = flipud(AstC(Ic).Cat);
    AstC(Ic).SortedBy    = [];
    AstC(Ic).SortedByCol = [];
end