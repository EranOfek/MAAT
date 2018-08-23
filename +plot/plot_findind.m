function [AllInd,Data,Select]=plot_findind(Data,varargin);
%---------------------------------------------------------------------------
% plot_findind function                                            plotting
% Description: Select data points on a plot using the mouse
%              Return the data points and indices.
% Input  : - Optional [X,Y] to plot.
%            If not given or empty (i.e., []), then get the data points
%            from the current axis.
%          * Arbitrary number of arguments to pass to the plot command.
% Output : - Indices of selected data points.
%          - [X, Y] of the input data points
%          - The mouse selected positions.
% Tested : Matlab 7.8
%     By : Eran O. Ofek               November 2009
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [AllInd,Data,Select]=plot_findind(rand(100,2));
%---------------------------------------------------------------------------
import Util.Geom.*

MarkPlotSymbol = 'rx';
MarkPlot       = 'y';


Def.Data = [];

if (nargin==0),
   Data = Def.Data;
end


if (isempty(Data)==1),
   % get data from current axis
   X = get(get(gca,'Children'),'XData');
   Y = get(get(gca,'Children'),'YData');
   Data = [X.',Y.'];
else
   X = Data(:,1);
   Y = Data(:,2);
   if (nargin>1),
      plot(X,Y,varargin{:});
   else
      plot(X,Y,'o');
   end
end

IsHold = ishold;
hold on;

fprintf('Select data points using left mouse button, abort using right click\n');

B = 1;
fprintf('Selected points:\n');
fprintf('N, I, X, Y, Dist:\n');
Counter = 0;
AllInd = zeros(0,1);
Select = zeros(0,2);
while (B==1),
   Counter = Counter + 1;
   [Xp,Yp,B]=ginput(1);

   % select nearest point
   Dist = plane_dist(Xp,Yp,X,Y);
   [MinDist,MinInd] = min(Dist);

   switch lower(MarkPlot)
    case 'y'
       plot(X(MinInd),Y(MinInd),MarkPlotSymbol);
    otherwise
       % do nothing
   end

   fprintf('%5d %9d   %f %f   %f\n',Counter,MinInd, X(MinInd),Y(MinInd),MinDist);
   AllInd = [AllInd; MinInd];
   Select = [Select; [Xp, Yp]];
   
end


switch IsHold
 case 0
    hold off
 case 1
    hold on
 otherwise
    error('Unknown IsHold option');
end
