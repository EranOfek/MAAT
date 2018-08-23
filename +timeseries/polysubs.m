function [y,a]=polysubs(x,n,c_x,c_y)
% Fit and subtract polynomials from a timeseries [T,Y].
% Package: timeseris
% Description: Subtract polynomial from a data set (no errors).
% Input  : - Data matrix.
%          - Degree of polynom to subtract.
%          - c_x, column of dependent variable, default is 1.
%          - c_y, column of independent variable, default is 2.
% Output : - Data matrix after the polynomial subtraction.
%            The other columns of original matrix are copied to y
%            only if c_x=1 and c_y=2.
%          - Polynomial coefecient of fit. 
% Tested : Matlab 4.2
%     By : Eran O. Ofek                       May 1994
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

if nargin==2
   c_x = 1;
   c_y = 2;
elseif nargin==3
   c_y = 2;
elseif nargin==4
   % do nothing
else
   error('Illegal number of input arguments');
end

wid = length(x(1,:));
a = polyfit(x(:,c_x),x(:,c_y),n);
temp = [x(:,c_x), x(:,c_y) - polyval(a,x(:,c_x))];
if wid>c_y
   y = [temp,x(:,c_y+1:wid)];
else
   y = temp;
end
