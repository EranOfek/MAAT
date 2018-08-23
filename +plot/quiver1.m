function Harrow=quiver1(XYZ,UVW,varargin)
%--------------------------------------------------------------------------
% quiver1 function                                                plotting
% Description: An improved quiver function. Allows to control the arrows 
%              properties and shape.
% Input  : - Two or three column matrix containing [X, Y] or [X, Y, Z]
%            coordinates of the arrow base point.
%          - Two or three column matrix containing [X, Y] or [X, Y, Z]
%            coordinates of the arrow length component in each axis.
%          * Arrow properties. Pairs of ...,keyword,value,... containing
%            one or serveral of the following keywords:
%            'Length'       - Length of arrow head, see arrow for details.
%            'BaseAngle'    - See arrow for details.
%            'TipAngle'     - See arrow for details.
%            'Width'        - Arrow width, see arrow for details.
%            'CrossDir'     - See arrow for details.
%            'NormalDir'    - See arrow for details.
%            'LineStyle'    - Arrow base line style {'-' | '--' | '-.' | ':'}
%            'FaceColor'    - Arrow head face color.
%            'EdgeColor'    - Arrow edge color.
% Output : - Vector of handles for each arrow.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Mar 2004
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

NextPlot = get(gcf,'NextPlot'); %--- Hold handling ---
hold on;
box on;

N1 = size(XYZ,2);
N2 = size(UVW,2);
if (N1==2 && N2==2),
   %--- 2d quiver ---
   X  = XYZ(:,1);
   Y  = XYZ(:,2);
   U  = UVW(:,1);
   V  = UVW(:,2);

   Xh = X + U;   % X head coordinates
   Yh = Y + V;   % Y head coordinates
   %--- find best axis ---
   MinX = min([X;Xh]);
   MaxX = max([X;Xh]);
   MinY = min([Y;Yh]);
   MaxY = max([Y;Yh]); 
   
   
   axis([MinX MaxX MinY MaxY]);
   
   N = length(X);
   
   Harrow = zeros(N,1);
   for I=1:1:N,
      axis(axis)
      H=arrow([X(I);Y(I)],[Xh(I);Yh(I)],varargin{:});
      Harrow(I) = H;
   end

elseif (N1==3 && N2==3),
   %--- 3d quiver ---
   X  = XYZ(:,1);
   Y  = XYZ(:,2);
   Z  = XYZ(:,3);
   U  = UVW(:,1);
   V  = UVW(:,2);
   W  = UVW(:,3);

   Xh = X + U;   % X head coordinates
   Yh = Y + V;   % Y head coordinates
   Zh = Z + W;   % Z head coordinates
   %--- find best axis ---
   MinX = min([X;Xh]);
   MaxX = max([X;Xh]);
   MinY = min([Y;Yh]);
   MaxY = max([Y;Yh]); 
   MinZ = min([Z;Zh]);
   MaxZ = max([Z;Zh]); 
   
   
   axis([MinX MaxX MinY MaxY MinZ MaxZ]);
   
   N = length(X);
   
   Harrow = zeros(N,1);
   for I=1:1:N,
      axis(axis)
      H=arrow([X(I);Y(I);Z(I)],[Xh(I);Yh(I);Zh(I)],varargin{:});
      Harrow(I) = H;
   end
else
   error('Illegal number of dimensions');
end
   
   
%--- return hold to original state ---
set(gcf,'NextPlot',NextPlot);
