function [Mag]=blackbody_mag_c(Temp,FilterSystem,FilterName,MagSys,Radius,Dist)
%--------------------------------------------------------------------------
% blackbody_mag_c function                                       AstroSpec
% Description: Calculate the magnitude, in a given filter, of a black-body,
%              given its temperature, radius and distance.
% Input  : - Vector of black body temperature [K].
%          - Filter system, see get_filter.m for details.
%            Alternatively, this can be a filter transmission curve
%            [Wave, Trans].
%          - Filter name, see get_filter.m for details.
%          - Magnitude system, 'AB' or 'Vega'.
%          - Black body radius [cm], default is 1cm.
%          - Distance [pc], default is 10pc.
% Output : - Black body observed magnitude.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Feb 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See Also: blackbody_mag.m, blackbody_bolmag.m
% Example : [Mag]=blackbody_mag_c([5770;8000],'Johnson','V','Vega',696000e5);
% Reliable: 1
%--------------------------------------------------------------------------
Pc  = 3.08567758e+18; % = get_constant('pc','cgs');

DefRadius = 1;
DefDist   = 10;
if (nargin==4),
   Radius = DefRadius;
   Dist   = DefDist;
elseif (nargin==5),
   Dist   = DefDist;
elseif (nargin==6),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (ischar(FilterSystem)),
    Filter = get_astfilter(FilterSystem,FilterName);
    WaveRange = (Filter(1).min_wl: (Filter(1).max_wl - Filter(1).min_wl)./100 :Filter(1).max_wl).';
    %WaveRange = [0.01;WaveRange;1e15];
else
    Norm = trapz(FilterSystem(:,1),FilterSystem(:,2));
    FilterSystem(:,2) = FilterSystem(:,2).*(1./Norm);
    WaveRange = (min(FilterSystem(:,1)):range(FilterSystem(:,1))./100:max(FilterSystem(:,1))).';
end
    
Nt  = length(Temp);

if (length(Radius)==1),
   Radius = Radius.*ones(Nt,1);
end
if (length(Dist)==1),
   Dist = Dist.*ones(Nt,1);
end

Mag    = zeros(Nt,1).*NaN;
for It=1:1:Nt,
   %[Il,If] = black_body(Temp(It),WaveRange);
   [Il] = AstSpec.blackbody(Temp(It),WaveRange);
   Spec    = [WaveRange, Il.Int .* 4.*pi.* Radius(It).^2 ./(4.*pi.*(Dist(It).*Pc).^2)];
   [Mag(It)] = synphot(Spec,FilterSystem,FilterName,MagSys);
end
