function [X,Y]=pr_sin(Theta,Phi,Center)
% Slant ortographic (SIN) projection
% Package: celestial.proj
% Description: Slant ortographic (SIN) projection
% Input  : - Longitude [rad].
%          - Latitude [rad].
%          - [Long, Lat] of projection center [rad].
% Output : - X
%          - Y
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [X,Y]=celestial.proj.pr_sin(1,1,[1 1])
% Reliable: 2
%--------------------------------------------------------------------------
RAD = 180./pi;

Theta_C = Center(1);
Phi_C   = Center(2);
Zeta = +cot(Theta_C).*sin(Phi_C);
Eta  = -cot(Theta_C).*cos(Phi_C);

X = +RAD.*(cos(Theta).*sin(Phi) + Zeta.*(1-sin(Theta)));
Y = -RAD.*(cos(Theta).*cos(Phi) - Eta .*(1-sin(Theta)));