function [SelColNames,ColInd,UseInd]=select_exist_colnames(AstC,ColNames)
%--------------------------------------------------------------------------
% select_exist_colnames function                             class/@AstCat
% Description: Given a cell array of column names, select the first column
%              name that appers in an AstCat object.
% Input  : - An AstCat object - The search will be done only on the first
%            element in the AstCat object.
%          - A column-ordered cell array of column names. Each row specify all the
%            possible names of a given column (the function will select the
%            first apperance that appears in the catalog).
%            The column name is selected seperatly from each column, so
%            multiple columns allow to deal with multiple names.
%            If vector of numbers than this will be treated as already
%            selected column indices.
% Output : - A cell array of the selected column names, one element per
%            column in the cell array of column names.
%          - A vector of column indices corresponding to the selected
%            column names.
%          - A vector of the line indices of the selected elements in the
%            cell array of column names.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [SelColNames,ColInd,UseInd]=select_exist_colnames(Sim,{'XWIN_IMAGE';'X'});
% Reliable: 
%--------------------------------------------------------------------------

if (ischar(ColNames))
    ColNames = {ColNames};
end

if (~iscell(ColNames))
    % Assume the user supplied the column indices in a two element matrix
    ColInd       = ColNames;
    SelColNames  = ind2colname(AstC,ColInd);
    UseInd       = 1;
else
    % The user supplied a cell array of column names
    % For each column in ColNames
    [~,Ncol] = size(ColNames);
    UseInd      = zeros(1,Ncol);
    SelColNames = cell(1,Ncol);
    for Icol=1:1:Ncol
        FoundInd    = find(is_colname(AstC(1),ColNames(:,Icol)));
        if (isempty(FoundInd))
            UseInd(Icol)      = NaN;
            SelColNames{Icol} = '';
        else
            UseInd(Icol)      = FoundInd(1);  % select the first
            SelColNames{Icol} = ColNames{UseInd,Icol};
        end
    end
    ColInd   = colname2ind(AstC(1),SelColNames);
end

    
