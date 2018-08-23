function [CurXY_2,CurXY_1,CurXY_0,varargout]=plot_rm(X,Y,Marker,CallFun,IndH,varargin);
%------------------------------------------------------------------------
% plot_rm function                                              plotting
% Description: Plot a 2-dimensional graph and interactively remove and
%              restore data points from the plot.
%              This function is being replaced by plot_int.m
% Input  : - Vector of X data to plot.
%          - Vector of Y data to plot.
%          - Plot marker, default is 'bo'. If empty, use default.
%          - Optional function to call where the 'f' option is
%            selected (see menu). The Function has the form:
%            [...] = Fun(X,Y); where X and Y are the current
%            X and Y which are not deleted and not marked by red cross.
%            Default is empty matrix (i.e., no function).
%          - The index of varargout of Fun in which a graphic handle
%            is stored. This graphic handle will be deleted before
%            the function is called again. Default is empty matrix
%            (i.e., do noting).
%          * Arbitrary number of input arguments to pass to function
%            to call.
% Output : - [X Y] list of points that were not deleted or marked with
%            red crosses.
%          - [X Y] list of points that were marked as red crosses.
%          - [X Y] list of points that were deleted.
%          * Arbitrary number of output arguments, which are the
%            output of the function to call (CallFun).
% Tested : Matlab 7.0
%     By : Eran O. Ofek                          May 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See also: plot_int.m
% Example: X = [1:1:20].';  Y=1+2.*X+0.5.*X.^2+randn(size(X)).*5;
%          [XY2,XY1,XY0,PolyPar]=plot_rm(X,Y,'bo','polyfit',[],2);
% Reliable: 2
%------------------------------------------------------------------------
import Util.array.*
import Util.Geom.*

DefMarker  = 'bo';
DefCallFun = [];
DefIndH    = [];

if (nargin==2),
   Marker  = DefMarker;
   CallFun = DefCallFun;
   IndH    = DefIndH;
elseif (nargin==3),
   CallFun = DefCallFun;
   IndH    = DefIndH;
elseif (nargin==4),
   IndH    = DefIndH;
end

if (isempty(Marker)==1),
   Marker   = DefMarker;
end

CrossedMarkerSize = 12;

GraphicH  = [];

H2 = plot(X,Y,Marker);
disp('Zoom is on')

HoldStatus = get(gcf,'NextPlot');
hold on;
zoom on;
plot_rm_menu;


X  = get(H2,'XData');
Y  = get(H2,'YData');
if (size(X,1)>size(X,2)),
   CurXY_2 = [X, Y];
else
   CurXY_2 = [X.', Y.'];
end
CurXY_1 = zeros(0,2);
CurXY_0 = zeros(0,2);
N  = length(X);
% Status: 2- in plot;  1- red cross; 0- deleted
Status = 2.*ones(N,1);   


Ans = input('Press a menu key to continue (or m for menu): ','s');
while (strcmpi(Ans,'q')==0),
   switch Ans
    case 'x'
       disp(' select a point using the mouse - and mark with red cross')
       [Xm,Ym] = ginput(1);
       Dist = plane_dist(Xm,Ym,CurXY_2(:,1),CurXY_2(:,2));
       [Min,MinInd] = min(Dist);

       CurXY_1 = [CurXY_1; CurXY_2(MinInd,:)];
       CurXY_2 = delete_ind(CurXY_2,round(MinInd));
       set(H2,'XData',CurXY_2(:,1),'YData',CurXY_2(:,2));

       % red cros
       H1 = plot(CurXY_1(:,1),CurXY_1(:,2),'rx','MarkerSize',CrossedMarkerSize);
    case 'X'
       disp(' select a point using the mouse - and delete')
       [Xm,Ym] = ginput(1);
       Dist = plane_dist(Xm,Ym,CurXY_2(:,1),CurXY_2(:,2));
       [Min,MinInd] = min(Dist);

       CurXY_0 = [CurXY_0; CurXY_2(MinInd,:)];
       CurXY_2 = delete_ind(CurXY_2,round(MinInd));
       set(H2,'XData',CurXY_2(:,1),'YData',CurXY_2(:,2));
    case 'r'
       disp(' select a rectangular region - and marke points with red cross')
       Rect = getrect;
       Xmin = Rect(1); Xmax = Rect(1)+Rect(3);
       Ymin = Rect(2); Ymax = Rect(2)+Rect(4);
       Id   = find(CurXY_2(:,1)>Xmin & CurXY_2(:,1)<Xmax & CurXY_2(:,2)>Ymin & CurXY_2(:,2)<Ymax);

       CurXY_1 = [CurXY_1; CurXY_2(Id,:)];
       CurXY_2 = delete_ind(CurXY_2,round(Id));
       set(H2,'XData',CurXY_2(:,1),'YData',CurXY_2(:,2));

       % red cros
       H1 = plot(CurXY_1(:,1),CurXY_1(:,2),'rx','MarkerSize',CrossedMarkerSize);
    case 'R'
       disp(' select a rectangular region - and marke points with red cross')
       Rect = getrect;
       Xmin = Rect(1); Xmax = Rect(1)+Rect(3);
       Ymin = Rect(2); Ymax = Rect(2)+Rect(4);
       Id   = find(CurXY_2(:,1)>Xmin & CurXY_2(:,1)<Xmax & CurXY_2(:,2)>Ymin & CurXY_2(:,2)<Ymax);

       CurXY_0 = [CurXY_0; CurXY_2(Id,:)];
       CurXY_2 = delete_ind(CurXY_2,round(Id));
       set(H2,'XData',CurXY_2(:,1),'YData',CurXY_2(:,2));

    case 'a'
       disp(' return the nearest red cross point')
       if (size(CurXY_1,1)==0),
          disp(' plot does not contain red croses - skip')
       else
          [Xm,Ym] = ginput(1);
          Dist = plane_dist(Xm,Ym,CurXY_1(:,1),CurXY_1(:,2));
          [Min,MinInd] = min(Dist);

          CurXY_2 = [CurXY_2; CurXY_1(MinInd,:)];
          CurXY_1 = delete_ind(CurXY_1,round(MinInd));
          set(H2,'XData',CurXY_2(:,1),'YData',CurXY_2(:,2));
          % delete the red cross from plot
          set(H1,'XData',CurXY_1(:,1),'YData',CurXY_1(:,2));
       end

    case 'b'
       disp(' return all red cross within rectangular area')
       if (size(CurXY_1,1)==0),
          disp(' plot does not contain red croses - skip')
       else
          Rect = getrect;
          Xmin = Rect(1); Xmax = Rect(1)+Rect(3);
          Ymin = Rect(2); Ymax = Rect(2)+Rect(4);
          Id   = find(CurXY_1(:,1)>Xmin & CurXY_1(:,1)<Xmax & CurXY_1(:,2)>Ymin & CurXY_1(:,2)<Ymax);

          CurXY_2 = [CurXY_2; CurXY_1(Id,:)];
          CurXY_1 = delete_ind(CurXY_1,round(Id));
          set(H2,'XData',CurXY_2(:,1),'YData',CurXY_2(:,2));
          % delete the red cross from plot
          set(H1,'XData',CurXY_1(:,1),'YData',CurXY_1(:,2));
       end

    case 'd'
       disp(' return all red cross')
       if (size(CurXY_1,1)==0),
          disp(' plot does not contain red croses - skip')
       else
          CurXY_2 = [CurXY_2; CurXY_1];
          CurXY_1 = zeros(0,2);
          set(H2,'XData',CurXY_2(:,1),'YData',CurXY_2(:,2));
          % delete the red cross from plot
          set(H1,'XData',CurXY_1(:,1),'YData',CurXY_1(:,2));
       end

    case 'c'
       disp(' return all the deleted and red cross points')
       if (size(CurXY_1,1)==0 & size(CurXY_0,1)==0),
          disp(' plot does not contain red croses or deleted points - skip')
       else
          CurXY_2 = [CurXY_2; CurXY_1; CurXY_0];
          CurXY_1 = zeros(0,2);
          CurXY_0 = zeros(0,2);
          set(H2,'XData',CurXY_2(:,1),'YData',CurXY_2(:,2));
          % delete the red cross from plot
          if (exist('H1')==1),
             set(H1,'XData',CurXY_1(:,1),'YData',CurXY_1(:,2));
          end
       end

    case 'A'
       disp(' Add artifial point')
       [ArtX, ArtY] = ginput(1);
       CurXY_2 = [CurXY_2;[ArtX, ArtY]];
       plot(ArtX,ArtY,Marker);


    case 'f'
       % deleted GraphicH (from last call of function)
       if (isempty(GraphicH)==1),
          % do nothing
       else
          delete(GraphicH);
       end

       % call function
       if (isempty(CallFun)==1),
          disp('No function to call')
       else
          N = nargout-3;
          varargout = cell(N,1);
          ParStr = '[';
          for I=1:1:N-1,
             ParStr = sprintf('%s varargout{%d},',ParStr,I);
          end
          ParStr = sprintf('%s varargout{%d}]',ParStr,N);
          eval(sprintf('%s = feval(CallFun,CurXY_2(:,1),CurXY_2(:,2),varargin{:});',ParStr));
       end

       if (isempty(IndH)==1),
          GraphicH = [];
       else
          GraphicH = varargout{IndH};
       end

    case 'm'
       disp('');
       disp('Display Menu:');
       plot_rm_menu;

    otherwise
       disp('');
       disp('Unknown option - options are:');
       plot_rm_menu;
   end
   Ans = input('Press a menu key to continue (or m for menu): ','s');

end


set(gcf,'NextPlot',HoldStatus);



function plot_rm_menu;
%---------------------
disp('----------')
disp('Menu:     ')
disp('    q - quit')
disp('    x - put a red cross on the nearest point')
disp('    X - delete the nearest point')
disp('    r - put a red cross on all the points inside a rectangular box')
disp('    R - delete all the points inside a rectangular box')
disp('    a - return the nearest point with a red cross')
disp('    b - return all the points with a red cross within a rectabgular region')
disp('    d - return all the red croses')
disp('    c - return all the deleted and red cross points')
disp('    A - Add artificial point')
disp('    f - call the function')
disp('    m - display this menu')
disp('----------')
