function [varargout]=colname2ind(AstC,ColName,RetNaN)
%--------------------------------------------------------------------------
% colname2ind function                                       class/@AstCat
% Description: Given AstCat object convert column name to
%              column index.
% Input  : - A single AstCat object.
%          - Colum name, or a cell array of column names.
%            If numeric then return the input as is.
%          - Flag indicating if to return NaN if colname doesn't exist.
%            Default is true.
% Output : * A vector of column indices.
%            If multiple output argument, then each output
%            argument cooresponds to one column in ColName.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ColInd=colname2ind(AstC,'XWIN_IMAGE')
%          ColInd=colname2ind(AstC,{'XWIN_IMAGE','YWIN_IMAGE'})
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<3),
    RetNaN = true;
end

if (numel(AstC)>1),
    error('colname2ind works on a single element AstCat');
end

if (~isnumeric(ColName)),
    % ColName is a string
    if (~iscell(ColName)),
        ColName = {ColName};
    end
    Ncn = numel(ColName);
    ColInd = zeros(Ncn,1);
    for Icn=1:1:Ncn,
        if (RetNaN),
            if (isfield(AstC.Col,ColName{Icn})),
                ColInd(Icn) = AstC.Col.(ColName{Icn});
            else
                ColInd(Icn) = NaN;
            end
        else
            ColInd(Icn) = AstC.Col.(ColName{Icn});
        end
    end
else
    % return as is
    ColInd = ColName;
end

if (nargout>1),
    % the user requested for the ColInd in
    % different parameters
    Ncn = numel(ColInd);
    CellColInd = num2cell(ColInd);
    [varargout{1:Ncn}] = deal(CellColInd{:});

else
    varargout{1} = ColInd;
end

