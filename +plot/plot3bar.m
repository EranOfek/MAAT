function plot3bar(x,y,z,symbol);
%-----------------------------------------------------------------------
% plot3bar function                                            plotting
% Description: plot 3-Dim graph with line connecting the points to
%              the z=0 surface.
% Input  : - vector of x data.
%          - vector of y data.
%          - vector of z data.
%          - the point symbol.
% Output : null
% Tested : Matlab 5.1
%     By : Eran O. Ofek           June 1997
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------------
N = length(x);
plot3(x,y,z,symbol);
hold on;

for i=1:1:N,
   plot3([x(i);x(i)],[y(i);y(i)],[0;z(i)]);
end
