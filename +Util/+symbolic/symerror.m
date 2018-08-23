function [ErrExp,ErrVar]=symerror(Expression,varargin)
% Calculate symbolic errors
% Package: Util.symbolic
% Description: Given a symbolic expression and the variables in the
%              expression, calculate the symbolic error function of the
%              expression, with respect to the variables. The output
%              error expression contains an error variable named
%              D_"original_var" for each variable in the input expression.
% Input  : - Symbolic expression.
%          * Arbitrary number of variables to differentiate by.
% Output : - Symbolic error expression.
%          - Cell array of new error variables.
% Tested : Matlab 6.5
%     By : Eran O. Ofek                    Aug 2004
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: symerror_calc.m
% Example: [ErrExp,ErrVar]=symerror('sqrt(x^2+y^2)','x','y')
%          % or
%          syms x y
%          [ErrExp,ErrVar]=symerror(sqrt(x^2+y^2),'x','y')
%          % or
%          syms x y
%          [ErrExp,ErrVar]=symerror(sqrt(x^2+y^2),x,y)
%          % Use: char(vectorize(ErrExp)) to convert the symbolic
%          % expression to a vectorized string that can be evaluated
%          % using eval. The function also returns a cell array of
%          % the new error variables (e.g.,  D_"original_var").
% Reliable: 2
%------------------------------------------------------------------------------
Nvar = length(varargin);
for I=1:1:Nvar
   Var         = varargin{I};
   if (ischar(Var))
      syms(Var)
   end
	ErrVar{I} = sym(sprintf('D_%s',char(Var)));
   SubDiff{I} = (diff(Expression,Var)*ErrVar{I})^2;
end

syms ErrExp;
ErrExp = 0;
for I=1:1:Nvar
   ErrExp = ErrExp + SubDiff{I};
end
ErrExp = sqrt(ErrExp); 
