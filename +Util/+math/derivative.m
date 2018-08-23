function Slope=derivative(X,Y)
% Numerical derivative of a row vector.
% Package: Util,math
% Description: Numerical derivative of a row vector.
% Input  : - X row vector.
%          - Y row vector.
% Output : - dY/dX numerical derivative.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Slope=Util.math.derivative(X,Y);
% Reliable: 2
%--------------------------------------------------------------------------


Xn = [X(1)-(X(2)-X(1)); X; X(end)+(X(end)-X(end-1))];
Yn = [NaN; Y; NaN];
YI = Util.interp.interp1_nan(Xn,Yn,'linear','extrap');

DX = diff(Xn);
DY = diff(YI);

Slope = (DY(1:end-1)./DX(1:end-1) + DY(2:end)./DX(2:end)).*0.5;