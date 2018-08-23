function L=lanczos(X,A)
%--------------------------------------------------------------------------
% lanczos function                                                 General
% Description: Calculate the Lanczos interpolation kernel function.
% Input  : - X
%          - Order. Default is 2.
% Output : - The Lanczos function evaluated at X.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: L=lanczos(X);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1),
    A =2;
end

L = sinc(X).*sinc(X./A);
L(abs(X)>A) = 0;

