function SubAstC=show(AstC,Cols,Rows,Disp)
% Display selected rows and columns of an AstCat object.
% Package: @AstCat
% Description: Display selected rows and columns of an AstCat object.
% Input  : - AstCat class object.
%          - Columns indices or a cell array of column names to display.
%            If empty, display all columns.
%          - Rows indices to display. This can be either:
%            1) a vector of rows indices;
%            2) a vector of logicals (true|false) indicating which
%               rows to display;
%            3) An empty matrix (display all rows).
%            4) A string with a number N - in this case the N last rows
%               will be displayed.
%            Default is empty.
%          - Display the AstCat {true|false}. Default is true.
% Output : - An AstCat class object with the selected columns.
% See also: AstCat/col_select.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstC=show(AstC,[],'10')
% Reliable: 2
%--------------------------------------------------------------------------

CatField   = 'Cat';

Def.Cols = [];
Def.Rows = [];
Def.Disp = true;
if (nargin==1)
    Cols = Def.Cols;
    Rows = Def.Rows;
    Disp = Def.Disp;
elseif (nargin==2)
    Rows = Def.Rows;
    Disp = Def.Disp;
elseif (nargin==3)
    Disp = Def.Disp;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments');
end


% display
Ncat    = numel(AstC);
if (SIM.issim(AstC))
    SubAstC = SIM(size(AstC));
else
    SubAstC = AstCat(size(AstC));
end

for Icat=1:1:Ncat
    
    % select
    if (isempty(Cols))
        ColInd = (1:1:size(AstC(Icat).(CatField),2));
    else
        ColInd = Cols;
    end
    
    SubAstC(Icat) = col_select(AstC(Icat),ColInd,Rows);
    
    if (Disp)
        if (~istable(SubAstC(Icat).(CatField)))
            DispSubAstC = astcat_array2table(SubAstC(Icat));
        end
    
        disp(DispSubAstC.(CatField));
    end
end

    

