function [AstC,Ind]=row_select(AstC,Ind)
% Select rows from a single element AstCat object.
% Package: @AstCat
% Description: Given a single AstCat object, return a new AstCat with
%              selected rows.
% Input  : - A single AstCat object.
%          - Indices or flags of rows to select.
%            Alternatively, a three column cell
%            array of {column_name, low_value, high_value}.
%            Only rows which column values are within these ranges will be
%            selected.
% Output : - An AstCat object with the selected rows.
%          - The vector of selected indices or logical flags.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstC=row_select(AstC,(6:10));
%          AstC=row_select(AstC,{'a',0,0.9;'b',0.5,1});
% Reliable: 2
%--------------------------------------------------------------------------

if (numel(AstC)>1)
    error('AstCat object must contain a single catalog');
end


if (iscell(Ind))
    % Ind is a cell array of constraints
    Ncrit = size(Ind,1);
    Flag  = true(sizecat(AstC),1);
    for Icrit=1:1:Ncrit
        Col  = col_get(AstC,Ind{Icrit,1});
        Flag = Flag & Col>Ind{Icrit,2} & Col<Ind{Icrit,3};
        
    end
    Ind = Flag;
end

% Ind is a vector of logicals or indices
AstC.Cat = AstC.Cat(Ind,:);
