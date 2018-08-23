function [AS1]=arith(AS1,AS2,Operator)
% Basic arithmetics on AstSpec class objects.
% Package: @AstSpec
% Description: Basic arithmetics on AstSpec class objects.
%              Binary operations on two AstSpec class objects, or
%              a AstSpec class object and a scalar.
% Input  : - AstSpec class object.
%          - AstSpec class object, a scalar or a vector.
%            If vector, then the vector size should be equal to the
%            size of the first argument.
%          - Operator. E.g, @plus.
% Output : - The result in an AstSpec class object.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: B = astspec_arith(A(1:2),A(1),@times);
% Reliable: 2
%--------------------------------------------------------------------------

%DefV. = 
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

%if (~isastspec(AS1)),
%    error('First argument must be of AstSpec class');
%end

N1 = numel(AS1);
N2 = numel(AS2);
if (~(N2==1 || N2==N1)),
    error('Size of second argument should be 1 or equal to that of the first argument');
end


if (~isastspec(AS2)),
    % assume AS2 is numeric
    AS2 = AS2.*ones(size(AS1));
    
    for I1=1:1:N1,
        I2 = I1;
        if (~isempty(AS1(I1).Int)),
            AS1(I1).Int  = Operator(AS1(I1).Int,AS2(I2));
        end
        if (~isempty(AS1(I1).Back)),
            AS1(I1).Back = Operator(AS1(I1).Back,AS2(I2));
        end
        if (~isempty(AS1(I1).Err)),
            if (any(strcmpi(func2str(Operator),{'times','rdivide'}))),
                AS1(I1).Err  = Operator(AS1(I1).Err,AS2(I2));
            else
                % do nothing to .Err
            end
        end
    end
else
    % AS2 is AstSpec
    N = max(N1,N2);
    for I=1:1:N,
        I1 = min(I,N1);
        I2 = min(I,N2);
        if (~isempty(AS1(I1).Int)),
            AS1(I1).Int = Operator(AS1(I1).Int,AS2(I2).Int);
        end
        if (~isempty(AS1(I1).Back)),
            if (~isempty(AS2(I2).Back)),
                AS1(I1).Back = Operator(AS1(I1).Back,AS2(I2).Back);
            else
                % Back of AS2 is not available - use AS1 back
                % do nothing
            end
        end
        if (~isempty(AS1(I1).Err)),
            if (~isempty(AS2(I2).Err)),
                if (any(strcmpi(func2str(Operator),{'times'}))),
                    [~,AS1(I1).Err] = times_err(AS1(I1).Int,AS1(I1).Err, AS2(I2).Int,AS2(I2).Err);
                elseif (any(strcmpi(func2str(Operator),{'rdivide'}))),
                    [~,AS1(I1).Err] = rdivide_err(AS1(I1).Int,AS1(I1).Err, AS2(I2).Int,AS2(I2).Err);
                elseif (any(strcmpi(func2str(Operator),{'plus','minus'}))),
                    AS1(I1).Err = sqrt(AS1(I1).Err.^2 + AS2(I2).Err.^2);
                else
                    % do nothing 
                end
            else
                % Err of AS2 is not available - use AS1 Err
                % do nothing
            end
        end
    end
    
end
        
        

        
    
    
    