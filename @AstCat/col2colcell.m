function [AstC,ColCell]=col2colcell(AstC)
%--------------------------------------------------------------------------
% col2colcell function                                       class/@AstCat
% Description: Repopulate the ColCell field based on the Col field.
% Input  : - AstCat object.
% Output : - AstCat object with the ColCell field repopulated.
%          - The ColCell field of the last element in the
%            AstCat array.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [AstC,ColCell]=col2colcell(AstC)
% Reliable: 2
%--------------------------------------------------------------------------




Nc = numel(AstC);
for Ic=1:1:Nc,
    Fields = fieldnames(AstC(Ic).Col);
    Vals   = struct2cell(AstC(Ic).Col);
    [~,SI] = sort(cell2mat(Vals));
    AstC(Ic).ColCell = Fields(SI).';
    ColCell = AstC(Ic).ColCell;
end

