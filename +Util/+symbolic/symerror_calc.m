function [FunVal,FunErr,ErrExp,ErrVar]=symerror_calc(Fun,varargin)
% Calculate and evaluate symbolic errors
% Package: Util.symbolic
% Description: Given a symbolic expression, the names of the variables
%              in the expression and the value of the variables and errors,
%              calculate the symbolic error function and evaluate it.
% Input  : - Symbolic function (e.g., 'sqrt(x^2+y^2)').
%          * Arbitrary triplets of arguments:
%            variable symbol followed by a matrix containing
%            the value and a second matrix containing the errors to
%            plug in each variable.
% Output : - Matrix of values of the evaluated expression.
%          - Matrix of errors of the evaluated expression.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Mar 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: symerror.m
% Example: syms x y
%          X = [1:1:10]'; ErrX = 0.1.*[1:1:10]';
%          Y = 6; ErrY = 0.2;
%          [Val,Err,ErrExp,ErrVar] = symerror_calc(sqrt(x^2+sin(y)), 'x',X,ErrX, 'y',Y,ErrY);
% Reliable: 2
%------------------------------------------------------------------------------

ColVal = 1;
ColErr = 2;

Narg = length(varargin);

% get the symbolic error function
VarInd = 0;
for I=1:3:Narg-2
   VarInd = VarInd + 1;
   Var{VarInd} = varargin{I};
   Val{VarInd} = varargin{I+ColVal};
   Err{VarInd} = varargin{I+ColErr};

   %eval(sprintf('%s=varargin{I+1};',varargin{I}));
end

[ErrExp,ErrVar] = Util.symbolic.symerror(Fun,Var{:});
ErrExpChar = vectorize(char(ErrExp));

VarInd = 0;
for I=1:3:Narg-2
   VarInd = VarInd + 1;
   
   eval(sprintf('%s = Val{VarInd};',char(Var{VarInd})));
   eval(sprintf('D_%s = Err{VarInd};',char(Var{VarInd})));
end

FunVal = eval(char(vectorize(Fun)));
FunErr = eval(ErrExpChar);


