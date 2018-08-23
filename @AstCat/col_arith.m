function [Val]=col_arith(AstC,Expression,OutType,ErrorH)
% Perform arithmetic operation on columns in an AstCat object.
% Package: @AstCat
% Description: Perform arithmetic operation on columns in an AstCat object.
% Input  : - AstCat object.
%          - An expression (string) or a cell array of expressions.
%            Each expression is evaluated and populate a column in the
%            output matrix. An expression contains column or columns names
%            and operations (see example).
%          - Output type:
%            'mat' - A matrix with column per expression. This output is
%                    valid only for single element AstCat object.
%            'astcat' - An AstCat object with element per catalog, and
%                    column per expression. Default.
%          - Error handling {true|false}. If false then will fail if
%            string can not be evaluated. If true, then will return
%            NaN if string can not be evaluated. Default is false.
% Output : - AstCat object or a matrix of evaluated expressions.
%            A column per expression.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Val = col_arith(AstC,'XWIN_IMAGE>5 & APER_MAG<18')
%          Val = col_arith(AstC,...
%               {'mod(XWIN_IMAGE,1)','XWIN_IMAGE+YWIN_IMAGE.^2','APER_MAG<18'});
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<3)
    OutType = 'astcat';
    ErrorH  = false;
elseif (nargin<4)
    ErrorH  = false;
else
    % do nothing
end

if (numel(AstC)>1 && ~strcmpi(OutType,'astcat'))
    error('OutType==mat works only for a single element AstCat');
end

if (~iscell(Expression))
    Expression = {Expression};
end

Nc = numel(AstC);
switch lower(OutType)
    case 'astcat'
        Val = AstCat(size(AstC));
    otherwise
        % do nothing
end
          
for Ic=1:1:Nc
    Nf = numel(AstC(Ic).ColCell); % number of fields
    if (istable(AstC(Ic).Cat))
        % table case
        for If=1:1:Nf
            % for each field
            % evaluate field into field name

            eval(sprintf(' %s = table2array(AstC(Ic).Cat(:,If));',AstC(Ic).ColCell{If}));
        end
    else
        % array case
        for If=1:1:Nf
            % for each field
            % evaluate field into field name

            eval(sprintf(' %s = AstC(Ic).Cat(:,If);',AstC(Ic).ColCell{If}));
        end
    end

    Nrow = size(AstC(Ic).Cat,1);
    % evaluate expressions
    Nexp = numel(Expression);   % number of expressions
    switch lower(OutType)
        case 'astcat'
            Val(Ic).Cat = zeros(Nrow,Nexp);
            Val(Ic).ColCell = cell(1,Nexp);
            
            if (ErrorH)
                % handle errors - replace with NaN
                for Iexp=1:1:Nexp
                    try
                        Val(Ic).Cat(:,Iexp) = eval(Expression{Iexp});
                    catch
                        Val(Ic).Cat(:,Iexp) = nan(Nrow,1);
                    end
                    Val(Ic).ColCell{Iexp} = sprintf('exp%d',Iexp);
                end
            else
                % do not handle errors
                for Iexp=1:1:Nexp
                    Val(Ic).Cat(:,Iexp) = eval(Expression{Iexp});
                    Val(Ic).ColCell{Iexp} = sprintf('exp%d',Iexp);
                end
            end
            Val(Ic) = colcell2col(Val(Ic));
        case 'mat'
            Val  = zeros(size(AstC(Ic).Cat,1),Nexp);
            if (ErrorH)
                % handle errors - replace with NaN
                for Iexp=1:1:Nexp
                    try
                        if (islogical(Expression{Iexp}))
                            Val(:,Iexp) = Expression{Iexp};
                        else
                            Val(:,Iexp) = eval(Expression{Iexp});
                        end
                    catch
                        Val(:,Iexp) = nan(Nrow,1);
                    end
                end
            else
                % do not handle errors    
                for Iexp=1:1:Nexp
                    if (islogical(Expression{Iexp}))
                        Val(:,Iexp) = Expression{Iexp};
                    else
                        Val(:,Iexp) = eval(Expression{Iexp});
                    end
                end
            end
        otherwise
            error('Unknown OutType option');
    end
end