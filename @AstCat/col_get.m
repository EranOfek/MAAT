function [Col,Units]=col_get(AstC,Col,Units,UnitsOut)
% Get column, in a vector format, from an AstCat object.
% Package: @AstCat
% Description: Get column, in a vector format, from an AstCat object.
% Input  : - An AstCat object with a single element.
%          - Either a string containing the column name,
%            a cell array of strings containing column names,
%            or a column index,
%            or a vector containing the column itself (will return this
%            column as is).
%          - Column units, to be returned if ColUnits field is not
%            populated. Default is empty.
%          - Output units. If provided then will convert column to this
%            units. If empty, thenn do nothing. Default is empty.
% Output : - Vector containing the requested column.
%          - String containing the column units. Empty if no units are
%            available.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Col,Units]=col_get(S,'XWIN_IMAGE');
%          [Col,Units]=col_get(S,1);
%          [Col,Units]=col_get(S,S.Cat(:,1));
%          [Col,Units]=col_get(S,'XWIN_IMAGE','pix');
% Reliable: 2
%--------------------------------------------------------------------------

Def.Units    = [];
Def.UnitsOut = [];
if (nargin==2)
    Units    = Def.Units;
    UnitsOut = Def.UnitsOut;
elseif (nargin==3)
    UnitsOut = Def.UnitsOut;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments');
end


if (numel(AstC)>1)
    error('col_get works on a single element AstCat object');
end

Nrow = size(AstC.Cat,1);

if (isnumeric(Col) && numel(Col)==Nrow)
    % assume Col is the column itself
    if (isempty(Units) && nargout>1)
        Units = [];
    end
    % do nothing
    % Col contains Col
    % Units contains Units
elseif (isnumeric(Col) && numel(Col)==1)
    % assume Col contains column index
    ColInd = Col;
    Col    = AstC.Cat(:,ColInd);
    if (isempty(Units))
        if (~isempty(AstC.ColUnits))
            Units = AstC.ColUnits{ColInd};
        else
            Units = [];
        end
    else
        % do nothing
        % Units contains Units
    end
elseif (ischar(Col) || iscellstr(Col))
    % Col contains column name
    ColInd = colname2ind(AstC,Col);
    if (isnan(ColInd))
        Col = NaN;
    else
        Col    = AstC.Cat(:,ColInd);
    end
    
    if (isempty(Units))
        if (~isempty(AstC.ColUnits))
            Units = AstC.ColUnits{ColInd};
        else
            Units = [];
        end
    else
        % do nothing
        % Units contains Units
    end
else
    error('Unknown Col input');
end

    

% convert to output units
if (~isempty(UnitsOut) && ~isempty(Units))
    Col = Col.*convert.units(Units,UnitsOut);
end
    
    
    
    