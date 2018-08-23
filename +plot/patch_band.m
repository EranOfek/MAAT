function H=patch_band(X,Y,WidthY,Color)
%--------------------------------------------------------------------------
% patch_band function                                             plotting
% Description: Given X and Y positions describing a line, plot a band
%              (i.e., using the patch command), around the line, where
%              the patch is defined to have a given height below and 
%              above the line.
% Input  : - X
%          - Y
%          - Semi width of patch in Y direction
%          - Color, default is [0.8 0.8 0.8].
% Output : - Patch handle
% Tested : Matlab 7.8
%     By : Eran O. Ofek                    Nov 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=[0:1:10]'; Y=X.^2; WidthY = 1;  H=patch_band(X,Y,WidthY);
% Reliable: 2
%--------------------------------------------------------------------------
Def.Color = [0.8 0.8 0.8];

if (nargin==3),
   Color = Def.Color;
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (length(WidthY)==1),
   WidthY = WidthY.*ones(size(X));
end

AllX = [X; flipud(X)];
AllY = [Y-WidthY; flipud(Y+WidthY)];

H = patch(AllX,AllY,Color);
set(H,'EdgeColor',Color);
