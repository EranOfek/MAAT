function Points=plot_select_points(X,Y,NearAlgo)
%--------------------------------------------------------------------------
% plot_select_points function                                     plotting
% Description: Given an open plot, and X/Y vectors, let the use select
%              using the left-mouse-click arbitrary number of points.
%              For each points return the clicked position, and the
%              nearest in the [X,Y] vectors.
%              Right-mouse click to terminate.
%              Use middle click to remove a point.
% Input  : - X vector.
%          - Y vector.
%          - Nearest point algorithm:
%            'X' - nearest in x position (default).
%            'D' - nearest in distance.
% Output : - A structure array of the selected points.
%            The following fields are available:
%            .ClickX  - The user X click position.
%            .ClickY  - The user Y click position.
%            .X       - X position of nearest point.
%            .Y       - Y position of nearest point.
%            .I       - Index of nearest point.
% Tested : Matlab R2011b
%     By : Eran Ofek / Ofer Yaron          Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=(1:1:100)'; Y=sin(2.*pi.*X./30); plot(X,Y)
%          Points=plot_select_points(X,Y);
% Reliable: 2
%--------------------------------------------------------------------------
import Util.Geom.*

PlotSymbol = 'bx';

if (nargin==2),
    NearAlgo = 'X';
end


Cont = true;
Ind = 0;
H = [];
while Cont
   Res = plot_int1(gcf,'mouse');
   
   
   switch Res.MB
       case 1
           
           Ind = Ind + 1;
           Points(Ind).ClickX = Res.Pos(1);
           Points(Ind).ClickY = Res.Pos(2);

           switch lower(NearAlgo)
               case 'd'
                   D = plane_dist(Res.Pos(1),Res.Pos(2),X,Y);
                   [~,MinI] = min(D);
               case 'x'
                   D = Res.Pos(1) - X;
                   [~,MinI] = min(abs(D));
               otherwise
                   error('Unknown NearAlgo option');
           end

           Points(Ind).X = X(MinI);
           Points(Ind).Y = Y(MinI);
           Points(Ind).I = MinI;
           Points(Ind).Status = 1;
           
           % plot point
           hold on;
           H(Ind) = plot(Points(Ind).X,Points(Ind).Y,PlotSymbol);

       case 3
           Cont = false;
           
       case 2
           % remove point
           D = plane_dist([Points.X]',[Points.Y]',Res.Pos(1),Res.Pos(2));
           [~,MinI] = min(D);
           Points(MinI).Status = 0;
           delete(H(MinI));
           
       otherwise
           % do nothing
   end
end

Points = Points([Points.Status]==1);
