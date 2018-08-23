function H=plot_ellipse(Position,Axes,Eccen,PA,Color,Width,CosDec,FaceColor)
%--------------------------------------------------------------------------
% plot_ellipse function                                           plotting
% Description: Plot an ellipse with given properties.
% Input  : - Position [X,Y] of ellipse center within current axis.
%          - [Semi Major axis, Semi Minor axis]
%          - Eccentricity. If given (e.g., not empty []) then Minor axis
%            is calculated from major axis and eccentricity.
%          - PA [rad], measured from X axis counterclocwise.
%          - Color, default is [0 0 0] (e.g., 'k', black).
%          - Line Width, default is 1.
%          - if cos(Dec) is given then stretch X axis by 1/cos(Dec),
%            default is 1.
%          - Ellipse face color, default is 'none'.
% Output : - Ellipse handle
% Plot   : - An ellipse plot
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Mar 2004
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: H=plot_ellipse([1,1],[1 0.5],[],1,'r',2);
% Reliable: 2
%--------------------------------------------------------------------------
RAD = 180./pi;

DefColor     = [0 0 0];
DefWidth     = 1;
DefCosDec    = 1;
DefFaceColor = 'none';

if (nargin==4),
   Color     = DefColor;
   Width     = DefWidth;
   CosDec    = DefCosDec;
   FaceColor = DefFaceColor;
elseif (nargin==5),
   Width     = DefWidth;
   CosDec    = DefCosDec;
   FaceColor = DefFaceColor;
elseif (nargin==6),
   CosDec    = DefCosDec;
   FaceColor = DefFaceColor;
elseif (nargin==7),
   FaceColor = DefFaceColor;
elseif (nargin==8),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(Eccen)==1),
   MajorAxis = Axes(1);
   MinorAxis = Axes(2);
else
   MajorAxis = Axes(1);
   MinorAxis = sqrt(MajorAxis(1).^2.*(1 - Eccen.^2));
end


Theta = (0:5:360)'./RAD;

X     = MajorAxis.*cos(Theta);
Y     = MinorAxis.*sin(Theta);

NX    = Position(1) + (cos(PA).*X - sin(PA).*Y)./CosDec;
NY    = Position(2) + sin(PA).*X + cos(PA).*Y;


%H     = plot(NX,NY);
H      = patch(NX,NY,'k');
set(H,'EdgeColor',Color,'FaceColor',FaceColor,'LineWidth',Width);

