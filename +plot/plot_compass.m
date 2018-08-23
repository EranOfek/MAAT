function [Hcompass,Htext]=plot_compass(Pos,Length,Color,Rotation,Dir);
%-----------------------------------------------------------------------
% plot_compass function                                        plotting
% Description: Plon a compass on a plot or image.
% Input  : - Compass position [X, Y].
%          - Compass arms length [pixels].
%          - Color, default is 'k';
%          - Angle of compass north in deg.
%            Default is 0. (north is up)
%          - Compass north east direction:
%            -1 : north-east direction is clockwise.
%            +1 : north-east direction is counter-clockwise.
%            Default is +1.
% Output : - Handle for the compass line.
%          - Handles for the text.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%-----------------------------------------------------------------------
RAD        = 180./pi;
DistFactor = 1.15;

if (nargin==2),
   Color     = 'k';
   Rotation  = 0;
   Dir       = +1;
elseif (nargin==3),
   Rotation  = 0;
   Dir       = +1;
elseif (nargin==4),
   Dir       = +1;
elseif (nargin==5),
   % do nothing
else
   error('Illegal Number of input arguments');
end

NextPlot = get(gca,'NextPlot');
hold on;

XLim   = get(gca,'XLim');
YLim   = get(gca,'YLim');
Xdiff  = abs(diff(XLim));
Ydiff  = abs(diff(YLim));


CompassX = Pos(1) + Length.*Dir.*[+sin(Dir.*Rotation./RAD);0;-cos(Dir.*Rotation./RAD)];
CompassY = Pos(2) + Length.*[+cos(Dir.*Rotation./RAD);0;+sin(Dir.*Rotation./RAD)];
TextX    = Pos(1) + DistFactor.*Length.*Dir.*[+sin(Dir.*Rotation./RAD);-cos(Dir.*Rotation./RAD)];
TextY    = Pos(2) + DistFactor.*Length.*[+cos(Dir.*Rotation./RAD);+sin(Dir.*Rotation./RAD)];

Hcompass = plot(CompassX,CompassY);
set(Hcompass,'Color',Color);


%--- plot text ---
Htext     = zeros(2,1);
Htext(1)  = text(TextX(1),TextY(1),'N');
Htext(2)  = text(TextX(2),TextY(2),'E');
set(Htext(1),'HorizontalAlignment','center','VerticalAlignment','middle')
set(Htext(2),'HorizontalAlignment','center','VerticalAlignment','middle')

set(gca,'XLim',XLim);
set(gca,'YLim',YLim);

set(gca,'NextPlot',NextPlot);
