function [Res]=fit_parabola2d_stationary(X,Y,Z,varargin)
% SHORT DESCRIPTION HERE
% Package: Util
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res]=Util.fit.fit_parabola2d_stationary
% Reliable: 
%--------------------------------------------------------------------------

if (nargin==0)
    % simulation / testing mode
   
    [X,Y] = meshgrid((1:1:100),(1:1:110));
    Z     = 1 + 0.3.*X + 0.2.*X.^2 - 0.5.*Y + 0.03.*Y.^2 + 0.02.*X.*Y; 
    Z     = 1 + (X-50).^2 + (Y-63).^2 + 0.02.*X.*Y;
    
    surface(X,Y,Z)
end

N = numel(Z);
H = [ones(N,1), X(:), X(:).^2, Y(:), Y(:).^2, X(:).*Y(:)];
Par = H\Z(:);

a0 = Par(1);
a1 = Par(2);
a2 = Par(3);
a3 = Par(4);
a4 = Par(5);
a5 = Par(6);


y_stat = (a3 - 2.*a2.*a5 - a1*a5)./(a5.^2 - 2.*a4);
x_stat = (-a3 - 2.*a4.*y_stat)./a5;

Type = 4.*a2.*a4 - a5.^2;

Res.Par = Par;
Res.x_stat = x_stat;
Res.y_stat = y_stat;
Res.Type   = Type;