function [Y,H]=lin_fun(Fun,Par,X,MapX)
% Evaluate a cell array of functions, linear in the free parameters.
% Package: Util.fit
% Description: Evaluate a cell array of functions, which are linear in the
%              free parameters:
%              i.e. Y=A(1)*Fun{1} + A(2)*Fun{2} + ...
% Input  : - Cell array of strings. Each cell containining a function
%            which multiply one of the coefficients.
%          - Column vector of coef.
%          - Vector of X (independent variable).
%            or matrix of independent variables if several of them
%            are needed (i.e. multi-dimensional fit; e.g. Y=f[X;Y;Z]).
%          - Cell vector in which the i-th element containing the column
%            in the X matrix which the i-th function is refering to.
%            if the i-th function has more than one dependent variable,
%            than the i-th cell should contain the column index (in the
%            X matrix) for each dependent variable,
%            by order of thie appearance.
%            If cell is empty, than function is a constant.
%            Default is assuming one-dim fit MapX{*}=[1];
% Output : - Value of function for each independent variable, X.
%          - Design matrix [Fun{1}, Fun{2}, ...];
% Tested : Matlab 7.0
%     By : Eran O. Ofek                July 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: % The function: Y=A + B*X^2 + C*sin(Y) + D*X*Y
%          N = 100;
%          X = [rand(N,1).*100, rand(N,1).*100];
%          Fun{1}='1'; Fun{2}='X.^2'; Fun{3}='sin(Y)'; Fun{4}='X.*Y';
%          MapX{1}=[]; MapX{2}=1; MapX{3}=2; MapX{4}=[1 2];
%          ParIn = [1; 0.1; 10; 0.02];
%          Y = lin_fun(Fun,ParIn,X,MapX);
% Reliable: 2
%------------------------------------------------------------------------------

Nfun = length(Fun);
Ny   = size(X,1);

if (nargin==3)
   for I=1:1:Ny
      MapX{I} = [1];
   end
elseif (nargin==4)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(MapX)==1)
  for I=1:1:Ny
	    MapX{I} = [1];
   end
end


H = zeros(Ny,Nfun);

% build the design matrix
for Ifun=1:1:Nfun
   InlineFun = inline(Fun{Ifun});

   Str = 'feval(InlineFun';

   if (isempty(MapX{Ifun})==1)
      Str = sprintf('%s,1',Str);
   else
      for Ivar=1:1:length(MapX{Ifun})
	 Str = sprintf('%s,X(:,MapX{Ifun}(:,%d))',Str,Ivar);
      end
   end
   Str = [Str,');'];

   H(:,Ifun) = eval(Str);
end

Y = sum(H*Par,2);
