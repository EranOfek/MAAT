function [Xc,Yc]=tri_equidist_center(x,y)
% Poistion of circumscribed circle.
% Package: Util.Geom
% Description: Find the position of a  point found in equal distances
%              from the verteces of a planer triangle.
% Input  : - Matrix of 3XN of triangles X coordinate verteces. Triangle per
%            line. 
%          - Matrix of 3XN of triangles Y coordinate verteces. Triangle per
%            line.
% Output : - Vector of triangle equi-distance X position from all verteces.
%          - Vector of triangle equi-distance Y position from all verteces.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Feb 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Xc,Yc]=tri_equidist_center([0;1;2],[1;0;4])
% Reliable: 2
%--------------------------------------------------------------------------

% solution is given by:
% syms x1 x2 x3 y1 y2 y3 xc yc
% A=solve('(xc-x1)^2+(yc-y1)^2=(xc-x2)^2+(yc-y2)^2','(xc-x2)^2+(yc-y2)^2=(xc-x3)^2+(yc-y3)^2',xc,yc)

Xc = (x(:,1).^2.*y(:,2) - x(:,2).^2.*y(:,1) - x(:,1).^2.*y(:,3) + x(:,3).^2.*y(:,1) + x(:,2).^2.*y(:,3) - x(:,3).^2.*y(:,2) - ...
    y(:,1).*y(:,2).^2 + y(:,1).^2.*y(:,2) + y(:,1).*y(:,3).^2 - y(:,1).^2.*y(:,3) - y(:,2).*y(:,3).^2 + y(:,2).^2.*y(:,3))./ ...
    (2.*(x(:,1).*y(:,2) - x(:,2).*y(:,1) - x(:,1).*y(:,3) + x(:,3).*y(:,1) + x(:,2).*y(:,3) - x(:,3).*y(:,2)));

Yc = (x(:,1).*x(:,2).^2 - x(:,1).^2.*x(:,2) - x(:,1).*x(:,3).^2 + x(:,1).^2.*x(:,3) + x(:,2).*x(:,3).^2 - x(:,2).^2.*x(:,3) + ...
    x(:,1).*y(:,2).^2 - x(:,2).*y(:,1).^2 - x(:,1).*y(:,3).^2 + x(:,3).*y(:,1).^2 + x(:,2).*y(:,3).^2 - x(:,3).*y(:,2).^2)./ ...
    (2.*(x(:,1).*y(:,2) - x(:,2).*y(:,1) - x(:,1).*y(:,3) + x(:,3).*y(:,1) + x(:,2).*y(:,3) - x(:,3).*y(:,2)));

