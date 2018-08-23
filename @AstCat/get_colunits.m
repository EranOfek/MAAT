function [ColUnits]=get_colunits(AstC,Col)
% Get the units for a specific column name or index in an AstCat object.
% Package: @AstCat
% Description: Get the units for a specific column name or index in an
%              AstCat object.
% Input  : - An AstCat object.
%          - Column name, cell array of column names, or column indices.
% Output : - A cell array of column units.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [ColUnits]=get_colunits(AstC,'XWIN_IMAGE');
% Reliable: 
%--------------------------------------------------------------------------

ColUnitsField   = 'ColUnits';

if (numel(AstC)>1),
    error('get_colunits except only a single element AstCat');
end


if (isempty(AstC.(ColUnitsField))),
    ColUnits = [];
else
    if (iscell(Col)),
        ColInd = colname2ind(AstC,Col);
    else
        ColInd = Col;
    end
    ColUnits = AstC.(ColUnitsField)(ColInd);
end

    
    