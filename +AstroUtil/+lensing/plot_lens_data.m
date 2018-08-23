function plot_lens_data(Theta,Beta,ColorIndex);
%-------------------------------------------------------------------------
% plot_lens_data function                                           glens
% Description: Given Image position and corresponding
%                             source position plot the image and source
%                             position with connecting lines.
% Input  : - Image position [X, Y].
%          - Source position [X, Y].
%          - Numerical index for color 1..17
% Output : null
% Plot   : Source and images position.
% Tested : Matlab 6.5
%     By : Eran O. Ofek     January 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-------------------------------------------------------------------------
ColX = 1;
ColY = 2;
ColorScheme    = [0.0  0.0  1.0
                  0.0  0.5  1.0
                  0.0  1.0  1.0
                  0.5  0.0  1.0
                  1.0  0.0  1.0
                  1.0  0.0  0.5
                  1.0  0.0  0.0
                  1.0  0.5  0.0
                  1.0  1.0  0.0
                  0.5  1.0  0.0
                  0.5  1.0  0.5
                  0.5  1.0  1.0
                  1.0  0.5  1.0
                  1.0  1.0  0.5
                  0.5  0.5  1.0
                  0.5  1.0  0.5
                  1.0  0.5  0.5];

if (nargin==2),
%   ColorScheme = ColorSchemeDef;
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end

Nt = length(Theta);

plot(Theta(:,ColX),Theta(:,ColY),'o','Color',ColorScheme(ColorIndex,:));
hold on;
plot(Beta(:,ColX),Beta(:,ColY),'.','Color',ColorScheme(ColorIndex,:));
for I=1:1:size(Theta,1),
   plot([Theta(I,ColX);Beta(I,ColX)],[Theta(I,ColY);Beta(I,ColY)],'-','Color',ColorScheme(ColorIndex,:));
end
