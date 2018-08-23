function [Ho,Xo,Yo]=gno(Hgca,varargin)
% Get the nearest object handle
% Package: plot
% Description: Get the nearest object handle. The user click near a
%              a figure objects, and the function return their handles
%              and coordinates.
%              Left click for object selection.
%              Right click for exiting program.
% Input  : - Handle to axis from which to select objects. If empty,
%            use current axis. Default is empty matrix.
%          * If empty than the nearest selected object will not be marked.
%            Otherwise, this can be symbol and color properties that
%            will be transfered to the plot function
%            (e.g., 'rx','MarkerSize',12). Default is empty matrix.
% Output : - Row vector of all objects handle.
%          - Row vector of X coordinates of selected objects.
%          - Row vector of Y coordinates of selected objects.
% Tested : Matlab 5.2
%     By : Eran O. Ofek                    Oct 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Ho,Xo,Yo]=gno;
%          [Ho,Xo,Yo]=gno(gca);
%          [Ho,Xo,Yo]=gno([],'rx','MarkerSize',12);
% Reliable: 2
%--------------------------------------------------------------------------
import Util.Geom.*

if (nargin==0)
    Hgca = [];
end

if (isempty(Hgca))
    Hgca   = gca;
end
% get objects positions
Children = get(Hgca,'Children');
Nc       = length(Children);

Data     = struct('Type',cell(1,Nc),'Xpos',cell(1,Nc),'Ypos',cell(1,Nc),'H',cell(1,Nc));
for I=1:1:Nc
    ObjType = get(Children(I),'Type');
    switch lower(ObjType)
        case 'line'
            Data(I).Type = ObjType;
            Data(I).Xpos = get(Children(I),'XData');
            Data(I).Ypos = get(Children(I),'YData');
            Data(I).H    = Children(I).*ones(size(Data(I).Xpos));
        case 'text'
            Data(I).Type = ObjType;
            Pos = get(Children(I),'Position');
            Data(I).Xpos = Pos(1);
            Data(I).Ypos = Pos(2);
            Data(I).H    = Children(I).*ones(size(Data(I).Xpos));
        otherwise
            error('Unknown Type');
    end
end
Xpos = [Data.Xpos].';
Ypos = [Data.Ypos].';
Hpos = [Data.H].';

Index = 0;
Button = 1;
Ho = [];
Xo = [];
Yo = [];
while (Button~=3)
   [X,Y,Button] = ginput(1);
   if (Button==1)
      Index = Index + 1;
      % search nearest object.
      D = plane_dist(X,Y,Xpos,Ypos);
      [~,MinI] = min(D);
      Ho(Index) = Hpos(MinI);
      Xo(Index) = Xpos(MinI);
      Yo(Index) = Ypos(MinI);
      
      % over plot
      if (~isempty(varargin))
          hold on;
          plot(Xo(Index),Yo(Index),varargin{:});
          hold off;
      end
   end
end
      

