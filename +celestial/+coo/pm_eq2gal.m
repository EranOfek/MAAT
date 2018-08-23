function [mu_l,mu_b]=pm_eq2gal(alpha,delta,mu_alpha,mu_delta,varargin)
% SHORT DESCRIPTION HERE
% Package: celestial
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [mu_l,mu_b]=celestial.coo.pm_eq2gal(alpha,delta,mu_alpha,mu_delta)
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

alpha_G = 192.85948./pi;
delta_G = -27.12825./pi;


DefV.CooUnits             = 'rad';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

CooFactor = convert.angular(InPar.CooUnits,'rad');
alpha = alpha.*CooFactor;
delta = delta.*CooFactor;

C1   = sin(delta_G).*cos(delta) - cos(delta_G).*sin(delta).*cos(alpha - alpha_G);
C2   = cos(delta_G).*sin(alpha - alpha_G);
Cosb = sqrt(C1.^2 + C2.^2);

Nsrc = numel(alpha);
mu_l = zeros(Nsrc,1);
mu_b = zeros(Nsrc,1);
for Isrc=1:1:Nsrc
    RotM = [C1(Isrc) C2(Isrc); -C2(Isrc) C1(Isrc)]./Cosb(Isrc);

    mu_G = RotM*[mu_alpha(Isrc); mu_delta(Isrc)];

    mu_l(Isrc) = mu_G(1,:);
    mu_b(Isrc) = mu_G(2,:);
end
