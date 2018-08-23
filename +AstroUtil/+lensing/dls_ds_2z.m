function Z_s=dls_ds_2z(Dls_Ds,Z_l,CosmoPars)
% Given D_ls/D_s ratio, z_l and cosmological parameters, solve for z_s.
% Package: AstroUtil.lensing
% Description: Given D_ls/D_s ratio, z_l and cosmological
%                        parameters, solve for z_s.
% Input  : - Vector of D_ls/D_s ratio.
%          - Vector of z_l
%          - Matrix of cosmological parameters [OmegaM, OmegaL].
% Output : - Vector of z_s.
% Tested : Matlab 6.5
%     By : Eran O. Ofek                    Feb 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%----------------------------------------------------------------------------

import AstroUtil.cosmo.*

H0 = 70; % doesn't change anything.
Tol   = 1e-3;

N = length(Dls_Ds);

Z_s = zeros(N,1);
for I=1:1:N
   Range = [Z_l(I)+Tol, 10];

   OmegaM = CosmoPars(I,1);
   OmegaL = CosmoPars(I,2);

   DlsDsFun = inline('ad_dist([Z_l Z_s],H0,OmegaM,OmegaL)./ad_dist([Z_s],H0,OmegaM,OmegaL)','Z_s','Z_l','H0','OmegaM','OmegaL');

   Z_s(I) = binsear_f(DlsDsFun,Dls_Ds(I),Range,Tol,Z_l(I),H0,OmegaM,OmegaL);
end
