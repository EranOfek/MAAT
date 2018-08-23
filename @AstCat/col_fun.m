function AstC=col_fun(AstC,Fun,Cols,Rows,AddPar)
% Evaluate a function that return a vector on columns in an AstCat object.
% Package: @AstCat
% Description: Evaluate a function that return a vector, on each column
%              in an AstCat object.
% Input  : - AstCat class object.
%          - The function to evaluate (e.g., @sin).
%          - Columns indices or a cell array of column names to evaluate.
%            If empty, display all columns.
%          - Rows indices to evaluate. This can be either:
%            1) a vector of rows indices;
%            2) a vector of logicals (true|false) indicating which
%               rows to display;
%            3) An empty matrix (display all rows).
%            4) A string with a number N - in this case the N last rows
%               will be displayed.
%            Default is empty.
%          - A cell array of additional parameters to pass to the function.
%            Default is {}.
% Output : - The input AstCat object in which the function was operated
%            on the requested columns and rows.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstC=col_fun(AstC,@sin)
%          % add 10 to the the first column
%          AstC=col_fun(AstC,@plus,1,[],10)
% Reliable: 2
%--------------------------------------------------------------------------


CatField         = 'Cat';

Def.Cols   = [];
Def.Rows   = [];
Def.AddPar = {};
if (nargin==2),
    Cols   = Def.Cols;
    Rows   = Def.Rows;
    AddPar = Def.AddPar;
elseif (nargin==3),
    Rows   = Def.Rows;
    AddPar = Def.AddPar;
elseif (nargin==4),
    AddPar = Def.AddPar;
elseif (nargin==5),
    % do nothing
else
    error('Illegal number of input arguments');
end

% if (~isempty(Rows) || ~isempty(Cols)),
%     AstC = show(AstC,Cols,Rows,false);
% end

[Nrow,Ncol] = sizecat(AstC);
Nc       = numel(AstC);
Val      = zeros(Nc,max(Ncol));
for Ic=1:1:Nc,
    
    if (ischar(Rows)),
        Irow = (max(1,Nrow(Ic)-str2double(Rows)):Nrow(Ic)).';
    else
        Irow = Rows;
    end
    if (isempty(Rows)),
        Irow = (1:1:Nrow(Ic)).';
    end
    
    if (isempty(Cols)),
        Icol = (1:1:Ncol(Ic));
    else
        Icol = colname2ind(AstC(Ic),Cols);
    end
    
    if (istable(AstC(Ic).(CatField))),
        % Catalog is table
        AstC(Ic).(CatField)(Irow,Icol) = Fun(table2array(AstC(Ic).(CatField)(Irow,Icol) ),AddPar{:});
    else
        AstC(Ic).(CatField)(Irow,Icol) = Fun(AstC(Ic).(CatField)(Irow,Icol),AddPar{:});
        
    end
end
