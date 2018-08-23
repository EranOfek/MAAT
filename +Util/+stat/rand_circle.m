function [X,Y,R]=rand_circle(varargin)
%------------------------------------------------------------------------------
% rand_circle function                                               AstroStat
% Description: Generate random number equally distributed inside a unit
%              circle.
% Input  : * Number of random numbers to generate:
%            e.g. (5,1)
% Output : - X.
%          - Y.
%          - R (radius).
% Tested : MATLAB 5.1
%     By : Eran O. Ofek                    Feb 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [X,Y,R]=rand_circle(10000,1);
%          plot(X,Y,'.');
% Reliable: 2
%------------------------------------------------------------------------------

T = 2*pi*rand(varargin{:});
U = rand(varargin{:})+rand(varargin{:});
R = U;
R(U>1) = 2-U(U>1);
X = R.*cos(T);
Y = R.*sin(T);
