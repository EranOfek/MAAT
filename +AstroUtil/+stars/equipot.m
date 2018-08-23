function [x,y,q]=equipot(m1,m2,a,s,n,z)
% Calculate the gravitational potential of a binary star on a grid.
% Package: AstroUtil.stars
% Description: Calculate two body equipotanials map.
% Input  : - m1 : first star mass in solar mass.
%          - m2 : second star mass in solar mass.
%          - a : the separation between the two stars. in meters.
%          - s : scaling the result from min to max in units of
%            the staller separation.
%          - n : number of point in each axis. default is 15.
%          - z : the surface to work on. default is z=0 (orbital plane)  
% Output : - grid of x coordinate.
%          - grid of y coordinate.
%          - matrix of potential defined by the x/y grids.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    May 1994
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [x,y,q]=AstroUtil.stars.equipot(1,0.7,0.7.*1e9,3,50);
%          mesh(x,y,q);
%-------------------------------------------------------------------------

if nargin==4,
   z = 0;
   n = 15;
elseif nargin==5,
   z = 0;
elseif nargin>6,
   error('4, 5 or 6 args only');
elseif nargin<4,
   error('4, 5 or 6 args only');
end
G = 6.672e-11;
m1 = m1.*1.9891e30;
m2 = m2.*1.9891e30;
mu = m2./(m1 + m2);
q  = m2./m2;
om2= G.*(m1 + m2)./(a.*a.*a);
d  = s.*a./n;
x  = (-s.*a:d:2.*s.*a);
y  = (-s.*a:d:2.*s.*a);
[MatX,MatY]=meshgrid(x,y);
t1 = G.*m1./sqrt(MatX.^2 + MatY.^2 + z.*z);
t2 = G.*m2./sqrt((MatX - a).*(MatX - a) + MatY.^2 + z.*z);
t3 = om2.*((MatX - mu.*a).*(MatX - mu.*a) + MatY.^2);
q  = -t1 - t2 - 0.5.*t3;
      