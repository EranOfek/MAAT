function [Val]=col_range(AstC,ColSelectRange,UnifyFun,OutType)
% Check if values in column is in a specified range.
% Package: @AstCat
% Description: Return logical flag indicating if in each row, the columns 
%              of the catalog are in specific range.
% Input  : - AstCat object.
%          - Three columns cell array indicating the column name
%            and its lower and upper range.
%            For example: {'APASS_htm_r',14,16;'APASS_htm_g',17,18}
%            will select entries with APASS r-band mag in the range 14 to
%            16 and APASS g-band mag in the range 17 to 18.
%          - A function by which to unify the selecttion criteria.
%            Empty matrix will do nothing and return a column of boolean
%            flags per criteria. @all will operate like the and operator,
%            and @any will operate like the or operator.
%            Default is @all.
%          - Output type:
%            'mat' - A vector with a logical flag. This output is
%                    valid only for single element AstCat object.
%            'astcat' - An AstCat object with element per catalog, and
%                    a single logical flag column. Default.
% Output : - AstCat object or a matrix of evaluated expressions.
%            A column per expression.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: B=col_range(A,{'XWIN_IMAGE',10,1000;'YWIN_IMAGE',10,1000});
% Reliable: 2
%--------------------------------------------------------------------------
Dim = 2;

if (nargin==2),
    UnifyFun = @all;
    OutType = 'astcat';
elseif (nargin==3),
    OutType = 'astcat';
elseif (nargin==4),
    % do nothing
else
    error('Illegal number of input arguments');
end

if (numel(AstC)>1 && ~strcmpi(OutType,'astcat')),
    error('OutType==mat works only for a single element AstCat');
end


Nc = numel(AstC);
switch lower(OutType)
    case 'astcat'
        Val = AstCat(size(AstC));
    otherwise
        % do nothing
end

Ncol = size(ColSelectRange,1);

for Ic=1:1:Nc,
    ColInd = colname2ind(AstC(Ic),ColSelectRange(:,1)');
    Ncat   = size(AstC(Ic).Cat,1);
    Flag = false(Ncat,Ncol);
    for Icol=1:1:Ncol,
        Flag(:,Icol) = AstC(Ic).Cat(:,ColInd(Icol))>=ColSelectRange{Icol,2} & ...
                       AstC(Ic).Cat(:,ColInd(Icol))>=ColSelectRange{Icol,2};
    end
    
    
    % evaluate expressions
    switch lower(OutType)
        case 'astcat'
            Val(Ic).Cat = UnifyFun(Flag,Dim);
            Val(Ic).ColCell{1} = 'flag';
            Val(Ic) = colcell2col(Val(Ic));
        case 'mat'
            Val  = UnifyFun(Flag,Dim);
        otherwise
            error('Unknown OutType option');
    end
end
