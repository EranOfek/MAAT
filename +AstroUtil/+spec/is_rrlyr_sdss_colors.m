function Flag=is_rrlyr_sdss_colors(u,g,r)
% Is RR Lyr star candidate based on SDSS colors
% Package: AstroUtil.spec
% Description: Select RR Lyr star candidates based on their SDSS magnitudes
%              in the u, g, and r-bands.
%              Based on 
% Input  : - u mag
%          - g mag
%          - r mag
% Output : - Flag indicating of RR Lyr star:
%            0 - no; 1 - 6%; 2 - 28%; 3 - 61%
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstroUtil.spec.is_rrlyr_sdss_colors(21,20,20)
% Reliable: 2
%--------------------------------------------------------------------------


%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

D_ug = (u - g) + 0.67.*(g - r) - 1.07;
D_gr = 0.45.*(u - g) - (g - r) - 0.12;

Flag1 = D_ug>-0.05 & D_ug<0.35 & D_gr>0.06 & D_gr<0.55;
Flag2 = D_ug>0.15  & D_ug<0.35 & D_gr>0.06 & D_gr<0.55;
Flag3 = D_ug>0.23  & D_ug<0.35 & D_gr>0.06 & D_gr<0.55;

Flag = zeros(size(u));
Flag(Flag1) = 1;
Flag(Flag2) = 2;
Flag(Flag3) = 3;

