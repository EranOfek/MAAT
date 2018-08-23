function H=plotline(X,Y,Len,Angle,Color,varargin);
%------------------------------------------------------------------------------
% plotline function                                                   plotting
% Description: Plot a line given the position of its start point, length,
%              and angle as measured in respect to the x-axis.
% Input  : - X (start)
%          - Y (start)
%          - Line length.
%          - Angle [deg], default is 90 (i.e. vertical line).
%          - Color and line type string, see plot for more details.
%            Default is 'b-',
%          * Additional arguments to pass to plot.
% Output : - Object handle.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                       May 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 1
%------------------------------------------------------------------------------
RAD = 180./pi;

if (nargin==3),
   Angle = 90;
   Color = 'b-';
elseif (nargin==4),
   Color = 'b-';
else
   % do nothing
end

N = length(X);
for I=1:1:N,
   H(I) = plot([X(I),X(I)+Len.*cos(Angle./RAD)],[Y(I),Y(I)+Len.*sin(Angle./RAD)],Color,varargin{:});
end
