function Mag=accretion_disk_mag_c(DiskPars,FilterSystem,FilterName,MagSys,Dist);
%------------------------------------------------------------------------------
% accretion_disk_mag_c function                                      AstroSpec
% Description: Calculate the magnitude, in a given filter, of a
%              optically-thick thin accretion disk model.
% Input  : - Matrix of accretion disk parameters:
%            [M Mdot Rin Rout]. See accretion_disk.m for details.
%          - Filter system, see get_filter.m for details,
%          - Filter name, see get_filter.m for details,
%          - Magnitude system, 'AB' or 'Vega'.
%          - Distance [pc], default is 10pc.
% Output : - Black body observed magnitude.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                       Feb 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See Also: blackbody_mag.m, accretion_disk.m
% Mag=accretion_disk_mag_c([1 1e-8 5 100],'SDSS','u','AB',10);
% Reliable: 2
%------------------------------------------------------------------------------
Pc  = get_constant('pc','cgs');

DefDist   = 10;
if (nargin==4),
   Dist   = DefDist;
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end


Filter = get_filter(FilterSystem,FilterName);
WaveRange = [Filter.min_wl{1}: (Filter.max_wl{1} - Filter.min_wl{1})./100 :Filter.max_wl{1}].';
WaveRange = [0.01;WaveRange;1e15];

Nt  = size(DiskPars,1);
Mag = zeros(Nt,1).*NaN;
for It=1:1:Nt,
   [Spec,VecR,VecT]=accretion_disk(DiskPars(It,1),DiskPars(It,2),DiskPars(It,3),DiskPars(It,4),WaveRange);
   Spec(:,2)  = Spec(:,2).*1e-8./(4.*pi.*(Dist.*Pc).^2);
   Mag(It) = synphot(Spec,FilterSystem,FilterName,MagSys);
end
