function AngRad=ang_radius_from_temp(G,Teff,A_G,Family,Band,System)
% Calculate the angular size of a star given its mag, extinction and temp
% Package: +AstroUtil.stars
% Input  : - Apparent Magnitude
%          - Effective temperature
%          - Extinction in band. Default is 0.
%          - Filter family. Default is 'GAIA'.
%          - Filter name. Default is 'G'.
%          - Magnitude system. Default is 'Vega'.
% Output : - Angular radius of star [arcsec]
% Example: [Cat,ColCell]=catsHTM.cone_search('GAIADR2',0,0,2031.1); 
% AngRad=AstroUtil.stars.ang_radius_from_temp(Cat(:,16),Cat(:,24),Cat(:,27))

RAD = 180./pi;

if nargin<6
    System = 'Vega';
    if nargin<5
        Band = 'G';
        if nargin<4
            Family = 'GAIA';
            if nargin<3
                A_G = 0;
            end
        end
    end
end

Gc = G - A_G;
[Mag]=AstroUtil.spec.blackbody_mag_c(Teff,Family,Band,System,constant.SunR,10,0);
SunAngRad = constant.SunR./(constant.pc.*10).*RAD.*3600;  % arcsec
AngRad = SunAngRad.*10.^(-0.2.*(Gc-Mag));   % arcsec

