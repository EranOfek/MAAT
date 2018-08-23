function Ind=plot_find_ind(Data,varargin);
%---------------------------------------------------------------------------
% plot_find_ind function                                           plotting
% Description: Plot data points and let the user select nearest points
%              to mouse position. Return the indices of selected points.
% Input  : (see errorxy.m for input options)
% Output : - Indices of selected points.
% Tested : Matlab 7.8
%     By : Eran O. Ofek                          September 2009
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------
import Util.Geom.*

CrossedMarkerSize = 12;

errorxy(Data,varargin{:});

HoldStatus = get(gcf,'NextPlot');
hold on;
zoom on;
plot_rm_menu;

X = Data(:,1);
Y = Data(:,2);
Ind = [];
Ans = input('Press a menu key to continue (or m for menu): ','s');
while (strcmpi(Ans,'q')==0),
   switch lower(Ans)
    case 'x'
       disp(' select a point using the mouse - and mark with red cross')
       [Xm,Ym] = ginput(1);
       Dist = plane_dist(Xm,Ym,X,Y);
       [Min,MinInd] = min(Dist);
       
       % plot red cross
       plot(X(MinInd),Y(MinInd),'rx','MarkerSize',CrossedMarkerSize)
       Ind = [Ind; MinInd];
    case 'r'
       disp(' select a rectangular region - and marke points with red cross')
       Rect = getrect;
       Xmin = Rect(1); Xmax = Rect(1)+Rect(3);
       Ymin = Rect(2); Ymax = Rect(2)+Rect(4);
       MinInd = find(X>Xmin & X<Xmax & Y>Ymin & Y<Ymax);

       % plot red cross
       plot(X(MinInd),Y(MinInd),'rx','MarkerSize',CrossedMarkerSize)
       Ind = [Ind; MinInd];
    case 'm'
       plot_rm_menu;
    otherwise
       plot_rm_menu;
   end
   Ans = input('Press a menu key to continue (or m for menu): ','s');

end

function plot_rm_menu;
%---------------------
disp('----------')
disp('Menu:     ')
disp('    q - quit')
disp('    x - put a red cross on the nearest point')
disp('    r - put a red cross on all the points inside a rectangular box')
disp('    m - display this menu')
disp('----------')
