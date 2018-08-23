function Val=col_fun2scalar(AstC,Fun,Cols,Rows)
%--------------------------------------------------------------------------
% col_fun2scalar function                                    class/@AstCat
% Description: Evaluate a function that return a scalar, on each column
%              in an AstCat object.
% Input  : - AstCat class object.
%          - The function to evaluate (e.g., @mean).
%          - Columns indices or a cell array of column names to display.
%            If Inf, display all columns. Default is Inf.
%          - Rows indices to display. This can be either:
%            1) a vector of rows indices;
%            2) a vector of logicals (true|false) indicating which
%               rows to display;
%            3) An empty matrix (display all rows).
%            4) A string with a number N - in this case the N last rows
%               will be displayed.
%            Default is empty.
% Output : - A matrix in which each row represent the mean value in each
%            column of an AstCat element.
% See also: AstCat/disp.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: M=col_fun2scalar(AstC,@mean)
% Reliable: 2
%--------------------------------------------------------------------------


CatField         = 'Cat';


Def.Cols = Inf;
Def.Rows = [];
if (nargin==2),
    Cols = Def.Cols;
    Rows = Def.Rows;
elseif (nargin==3),
    Rows = Def.Rows;
elseif (nargin==4),
    % do nothing
else
    error('Illegal number of input arguments');
end

% if (~isempty(Rows) || ~isempty(Cols)),
%     AstC = show(AstC,Cols,Rows,false);
% end
AstC = col_select(AstC,Cols,Rows);

[~,Ncol] = sizecat(AstC);
Nc       = numel(AstC);
Val      = zeros(Nc,max(Ncol));
for Ic=1:1:Nc,
    
    if (istable(AstC(Ic).(CatField))),
        % Catalog is table
        Val(Ic,:) = Fun(table2array(AstC(Ic).(CatField)));
    else
        Val(Ic,:) = Fun(AstC(Ic).(CatField));
        
    end
end
