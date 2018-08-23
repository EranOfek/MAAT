function []=mildrel_synchrotron(varargin)
% SHORT DESCRIPTION HERE
% Package: AstroUtil.supernova
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV. = 
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


%%
gamma = @(beta) 1./sqrt(1-beta.^2);

beta = (0.1:0.01:0.9)';



gb   = gamma(beta).*beta;

Mgb  = 2.5e-6.*gb.^-7.8;
Mgb(gb<0.3) = 2.5e-6.*0.3.^-7.8 .* (gb(gb<0.3)./0.3).^-3.0;

n = 0.05;
Eps_B = 1./30;
Eps_e = 1./3;
gamma_max_gamma_min = 1e5;

nu = 1.4e9;

t = 100.*86400.*(Mgb./(gamma(beta).*4./3.*pi.*n)  .* ((1-beta)./beta).^3).^(1./3);  % [s]

Mt = Mgb./gamma(beta);


nu_c  = 1.6e17 .*((1-beta).^2./(beta.^3)) .*(30.*Eps_B).^(-3./2) .*(t./(100.*86400)).^-2;   % Hz
% for nu>nu_c
nuL_nuH = 4.*pi.*Eps_e.*gamma(beta).^2.*beta.^5.*n.*constant.mp.*constant.c.^5.*t.^2./(2.*log(gamma_max_gamma_min).*(1-beta).^3);  % [erg/s]
% nu<nu_c
nuL_nuL = 4.4e35.*(3.*Eps_e) .* gamma(beta).^2.*beta.^(13./2)./((1-beta).^4) .*(30.*Eps_B).^(3./4) .* (n./0.01).^(7./4) .* (t./(100.*86400)).^3 .* (nu./1e9).^(1./2);  % [erg/s]

nu_m = 3e7.*(3.*Eps_e).^2 .*(30.*Eps_B).^(1./2) .* (n./0.01).^(1./2) .*gamma(beta).^5.*beta.^4;  % [Hz]
nu_a = 4e7.*(30.*Eps_B).^(1./3).*(3.*Eps_e).^(1./3) .* gamma(beta).^2.*beta.^(5./3) ./( (1-beta).^(1./3) ) .* (n./0.01).^(2./3) .* (t./(100.*86400)).^(1./3); % [Hz]



%loglog(gb,Mgb)
loglog(t./86400,nuL_nuL./nu./(4.*pi.*(3.08e18.*40e6).^2)./1e-26)   % mJy



