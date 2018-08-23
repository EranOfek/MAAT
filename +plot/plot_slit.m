function H=plot_slit(Pos,Length,Width,Rotation,EdgeColor,FaceColor,OpenSlit);
%-----------------------------------------------------------------------
% plot_slit function                                           plotting
% Description: Plot a slit (box) on a plot or image.
% Input  : - Slit center position [X, Y].
%          - Slit length [pixels].
%          - Slit Width [pixels].
%          - Slit position angle [deg], default is 0.
%          - Edge Color, default is 'k';
%          - Face Color, default is 'w';
%          - Slit is open {'y' | 'n'}, default is 'n'.
% Output : - Handle for slit patch.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%-----------------------------------------------------------------------
RAD        = 180./pi;

if (nargin==3),
   Rotation  = 0;
   EdgeColor = 'k';
   FaceColor = 'w';
   OpenSlit  = 'n';
elseif (nargin==4),
   EdgeColor = 'k';
   FaceColor = 'w';
   OpenSlit  = 'n';
elseif (nargin==5),
   FaceColor = 'w';
   OpenSlit  = 'n';
elseif (nargin==6),
   OpenSlit  = 'n';
elseif (nargin==7),
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

CosRot = cos(Rotation./RAD);
SinRot = sin(Rotation./RAD);

Arm1Center   = Pos + 0.5.*Width.*[+CosRot, -SinRot];
Arm2Center   = Pos + 0.5.*Width.*[-CosRot, +SinRot];
Arm1X        = Arm1Center(1) + 0.5.*Length.*[+SinRot; -SinRot];
Arm1Y        = Arm1Center(2) + 0.5.*Length.*[+CosRot; -CosRot];
Arm2X        = Arm2Center(1) + 0.5.*Length.*[+SinRot; -SinRot];
Arm2Y        = Arm2Center(2) + 0.5.*Length.*[+CosRot; -CosRot];

Arm1Xextra   = Arm1Center(1) + 0.6.*Length.*[+SinRot; -SinRot];
Arm1Yextra   = Arm1Center(2) + 0.6.*Length.*[+CosRot; -CosRot];
Arm2Xextra   = Arm2Center(1) + 0.6.*Length.*[+SinRot; -SinRot];
Arm2Yextra   = Arm2Center(2) + 0.6.*Length.*[+CosRot; -CosRot];


switch OpenSlit
 case 'y'
    H=plot(Arm1Xextra,Arm1Yextra,':');
    set(H,'Color',EdgeColor);
    H=plot(Arm2Xextra,Arm2Yextra,':');
    set(H,'Color',EdgeColor);

    H=plot(Arm1X,Arm1Y,'-');
    set(H,'Color',EdgeColor);
    H=plot(Arm2X,Arm2Y,'-');
    set(H,'Color',EdgeColor);

 case 'n'

    X            = [Arm1X; flipud(Arm2X)];
    Y            = [Arm1Y; flipud(Arm2Y)];
    H            = patch(X,Y,EdgeColor);
    set(H,'FaceAlpha',0,'FaceColor',FaceColor,'EdgeColor',EdgeColor);
 otherwise
    error('Unknown Method Option');
end

set(gca,'XLim',XLim);
set(gca,'YLim',YLim);

set(gca,'NextPlot',NextPlot);
