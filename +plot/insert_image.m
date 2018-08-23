function H=insert_image(X,Y,Image,ImageSize,DispRange,Pos,AlphaImage)
% Insert an image to a matlab figure in a given position.
% Package: plot
% Description: Insert an image to a matlab figure in a given position.
% Input  : - X image position.
%          - Y image position.
%          - Image name or matrix.
%          - Image Siez [X Y].
%          - Display range [LOW HIGH], default is [0 1].
%          - Position scheme:
%            'cen'  : place the image center at the X/Y position (default).
%            'cor'  : place the image bottom left corner at the X/Y position.
%          - Optional alpah image, default is no alpha image.
% Output : - Handle for the axis containing the image.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%-------------------------------------------------------------------------

if (nargin==4)
   DispRange = [0 1];
   Pos       = 'cen';
   AlphaImage = NaN;
elseif (nargin==5)
   Pos       = 'cen';
   AlphaImage = NaN;
elseif (nargin==6)
   AlphaImage = NaN;
elseif (nargin==7)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (ischar(Image)==1)
   ImMat = imread(Image);
else
   ImMat = Image;
end

% get hold status
HoldStatus = get(gcf,'NextPlot');
hold on;

% get current axes
CurAxes    = get(gcf,'CurrentAxes');


XLim    = get(gca,'XLim');
YLim    = get(gca,'YLim');
AxesPos = get(gca,'Position');

% calculate absolute position of image
Xabs     = AxesPos(1)+AxesPos(3).*(X - XLim(1))./range(XLim);
Yabs     = AxesPos(2)+AxesPos(4).*(Y - YLim(1))./range(YLim);
Xabssize = AxesPos(3).*ImageSize(1)./range(XLim);
Yabssize = AxesPos(4).*ImageSize(2)./range(YLim);

switch Pos
 case 'cor'
    OffsetX = 0;
    OffsetY = 0;
 case 'cen'
    OffsetX = -0.5.*Xabssize;
    OffsetY = -0.5.*Yabssize;
 otherwise
    error('Unknown Pos Option');
end

% create image axes
H   = axes('Position',[Xabs+OffsetX, Yabs+OffsetY, Xabssize, Yabssize]);
Him = imshow(ImMat,DispRange);
if (isnan(AlphaImage)==1)
   % do nothing
else
   set(Him,'AlphaData',AlphaImage);
end



% set current axes to old axes
set(gcf,'CurrentAxes',CurAxes);

% set hold status to original
set(gcf,'NextPlot',HoldStatus);
