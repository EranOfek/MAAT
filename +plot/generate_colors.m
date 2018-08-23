function Colors=generate_colors(N,ColorMapType)
% Generate equally spaced N colors from a given color map.
% Package: plot
% Description: Generate equally spaced N colors from a given color map.
% Input  : - Number of requested colors.
%          - Color map. Default is 'jet'.
% Output : - A matrix in which each row is an RGB color triplet.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Colors=generate_colors(5);
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin==1)
    ColorMapType = 'jet';
end

CM = colormap(ColorMapType);

Size = size(CM,1);
Step = Size./(N+1);

Pos = floor(1:Step:Size-Step);

Colors = CM(Pos,:);





