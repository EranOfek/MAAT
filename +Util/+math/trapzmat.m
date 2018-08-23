function Trapz=trapzmat(X,Y,Dim)
% Trapezoidal numerical integration on columns or rows of matrices.
% Package: Util.math
% Description: Trapezoidal numerical integration on columns or rows of
%              matrices.
%              Contrary to trapz.m, the X input for this function can
%              be a matrix.
% Input  : - X matrix.
%          - Y matrix.
%          - Dimension along to preform the integration, default is 1.
% Output : - Vector of integrals.
% Tested : Matlab 7.6
%     By : Eran O. Ofek                    Jun 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------
Def.Dim = 1;

if (nargin==2)
   Dim   = Def.Dim;
elseif (nargin==3)
   % do nothing
else
   error('Illrgal number of input arguments');
end

switch Dim
 case 1
    % do nothing
 case 2
    X = X.';
    Y = Y.';
 otherwise
    error('Unknwon Dim option');
end

Trapz = sum(abs(diff(X,1,1)).*0.5.*( Y(1:end-1,:) + Y(2:end,:) ) );
