function ColName=ind2colname(AstC,ColInd)
% Given AstCat object convert column index to column name.
% Package: @AstCat
% Description: Given AstCat object convert column index to
%              column name.
% Input  : - A single AstCat object.
%          - Vector of colum indices.
%            If a string or a cell array of strings reten as
%            is.
% Output : - Cell array of column names.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ColName=colname2ind(AstC,[1 2])
% Reliable: 2
%--------------------------------------------------------------------------

   
if (numel(AstC)>1),
    error('ind2colname works on a single element AstCat');
end

if (isnumeric(ColInd)),
    Ncn = numel(ColInd);
    ColName = cell(1,Ncn);
    for Icn=1:1:Ncn,
        ColName{Icn} = AstC.ColCell{ColInd(Icn)};
    end
else
    % input is string - return as is.
    ColName = ColInd;
end
        