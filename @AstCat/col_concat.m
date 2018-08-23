function AstC=col_concat(AstC1,AstC2,Col1,Col2)
%--------------------------------------------------------------------------
% col_concat function                                        class/@AstCat
% Description: Concat the columns of two AstCat objects which have
%              the same number of rows.
% Input  : - An AstCat object.
%          - An AstCat object in which each catalog have the same number of
%            rows as in the first AstCat object.
%          - Vector of column indices, or a cell array of column names
%            from the first AstCat object to concat. If empty, then
%            use all columns. Default is empty.
%          - Vector of column indices, or a cell array of column names
%            from the second AstCat object to concat. If empty, then
%            use all columns. Default is empty.
% Output : - AstCat object with the concated columns.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstC=col_concat(AstC1,AstC2,[],{'XWIN_IMAGE','YWIN_IMAGE'})
% Reliable: 
%--------------------------------------------------------------------------


if (nargin==2)
    Col1 = [];
    Col2 = [];
elseif (nargin==3)
    Col2 = [];
else
    % do nothing
end

Nc1 = numel(AstC1);
Nc2 = numel(AstC2);

Nc   = max(Nc1,Nc2);
AstC = AstC1;

for Ic=1:1:Nc
    Ic1 = min(Nc1,Ic);
    Ic2 = min(Nc2,Ic);
    
    if (isempty(Col1))
        % select all columns
        Ncol = size(AstC1(Ic1).Cat,2);
        ColInd1 = (1:1:Ncol);
    else
        ColInd1 = colname2ind(AstC1(Ic1),Col1);
    end
    if (isempty(Col2))
        % select all columns
        Ncol = size(AstC2(Ic2).Cat,2);
        ColInd2 = (1:1:Ncol);
    else
        ColInd2 = colname2ind(AstC2(Ic2),Col2);
    end
    
    AstC(Ic).Cat      = [AstC1(Ic1).Cat(:,ColInd1), AstC2(Ic2).Cat(:,ColInd2)];
    AstC(Ic).ColCell  = [AstC1(Ic1).ColCell(ColInd1), AstC2(Ic2).ColCell(ColInd2)];
    AstC(Ic)          = colcell2col(AstC(Ic));
    if (~isempty(AstC1(Ic1).ColUnits) && ~isempty(AstC2(Ic2).ColUnits))
        AstC(Ic).ColUnits = [AstC1(Ic1).ColUnits(ColInd1), AstC2(Ic2).ColUnits(ColInd2)];
    end
end
    