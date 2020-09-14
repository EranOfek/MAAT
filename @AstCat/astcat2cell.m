function [CellT,ColCell,Col]=astcat2cell(AS)
% Convert an astcat object to a cell array of matrix/tables
% Package: @AstCat
% Description: Convert an astcat object to a cell array of matrix/tables.
% Input  : - An AstCat object
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - Cell array of matrix/tables.
%          - ColCell of the first AstCat element.
%          - Col structure.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


CatField = AstCat.CatField;

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

N = numel(AS);
CellT = cell(size(AS));
for I=1:1:N
    CellT{I} = AS(I).(CatField);
end

ColCell = AS(1).ColCell;
Col     = AS(1).Col;

