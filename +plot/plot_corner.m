function H=plot_corner(X,Y,DX,DY,LineType,varargin);
%----------------------------------------------------------------------------
% plot_corner function                                              plotting
% Description: Plot corners which lines are parallel to the axes.
% Input  : - Vector of X points of the corners.
%          - Vector of Y points of the corners.
%          - Vector of X length of the corners.
%          - Vector of Y length of the corners.
%          - Line type, default is 'k-'.
%          * Arbitrary of number of input arguments to be passed to the
%            plot command.
% Output : - Vector of handels for each corner lines.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                  November 2010
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%----------------------------------------------------------------------------

N = length(X);
for I=1:1:N,
   H(I) = plot([X(I);X(I);X(I)+DX(I)],...
               [Y(I)+DY(I);Y(I);Y(I)],...
               LineType,varargin{:});
end
