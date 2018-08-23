function []=synchrotron(varargin)
% SHORT DESCRIPTION HERE
% Package: AstroUtil.supernova
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV. = 
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

syms kB T_B nu c R ne P_S nu_s B gamma_e beta_e S_nu d sigmaT e me

T_B = S_nu*d^2 * (c/nu)^2 / (2*pi*kB*R^2);
P_S = 4/3 * sigmaT * c *beta_e^2 *gamma_e^2 * B^2/(8*pi)
nu_s = gamma_e^2 * e * B/(2*pi*me*c)

2*pi* kB *T_B * (nu/c)^2  / (1/3 * R * ne * P_S/nu_s * (nu/nu_s)^(1/3))

Dist = 85e6.*constant.pc;
B    = 1;
S_nu = 20.*1e-26;
gamma_e = 1;

beta_e  = 0.5;

beta_e  = (0.03:0.01:0.99)';
gamma_e = 1./(sqrt(1-beta_e.^2));

nu      = 1.4e9;
Radius = 1e16;

nu_s    = gamma_e.^2 .*constant.e .*B./(2.*pi.*constant.me.*constant.c);
U_B     = B.^2./(8.*pi);
P_S     = 4./3.*constant.sigmaT.*constant.c.*beta_e.^2.*gamma_e.^2.*U_B;
T_B     = S_nu.*Dist.^2 .*(constant.c/nu).^2 ./(2.*pi.*constant.kB.*Radius.^2);
t_syn   = gamma_e.*constant.me.*constant.c.^2./P_S;
P_numax = P_S./nu_s;

%(9.*S_nu.*Dist.^2.*e)./(B.*ne.*Radius.^3.*beta_e.^2.*c.^2.*me.*sigmaT.*((2.*c.*me.*nu.*pi)./(B.*e.*gamma_e.^2)).^(1./3))
1./((9.*S_nu.*Dist.^2.*constant.e)./(B.*( 3.*Mass./(4.*pi.*constant.mp  )   ).*beta_e.^2.*constant.c.^2.*constant.me.*constant.sigmaT.*((2.*constant.c.*constant.me.*nu.*pi)./(B.*constant.e.*gamma_e.^2)).^(1./3)))



2.*pi.*constant.kB.*T_B.*(nu./constant.c).^2

Radius.*