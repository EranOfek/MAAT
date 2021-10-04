function [RA, Dec, RAantiSun, DecantiSun] = earthShadowCoo(JD, Dist, Args)
    % Calculate the J2000.0 equatorial coordinates of the Earth shadow at a given height
    % Input  : - JD (UT1 time scale).
    %          - Topocentric distance to point in shadow for which to
    %            calculate the position. Default is 42164 km.
    %          * ...,key,val,...
    %            'DistUnits' - Default is 'km'.
    %            'GeoPos' - [Lon, Lat, Height] must be in [rad, rad, m].
    %                   Default is [35 32 0]./(180./pi);  % [rad rad m]
    %            'RefEllipsoid' - Default is 'WGS84'.
    %            'OutUnitsDeg' - Output is in degrees. Default is true.
    % Output : - J2000.0 RA of shadow point.
    %          - J2000.0 Dec of shadow point.
    %          - J2000.0 RA of anti Sun direction.
    %          - J2000.0 Dec of anti Sun direction.
    % Author : Eran Ofek (Oct 2021)
    % Example: JD = 2451545 + (0:0.1:365)';
    %          [RA, Dec, RAas, Decas] = celestial.SolarSys.earthShadowCoo(JD, 'OutUnitsDeg',false);
    %          plot(JD, RAD.*celestial.coo.sphere_dist_fast(RA,Dec,RAas,Decas))
   
    arguments
        JD
        Dist                      = 42164;
        Args.DistUnits            = 'km';
        Args.GeoPos               = [35 32 0]./(180./pi);  % [rad rad m]
        Args.RefEllipsoid         = 'WGS84';
        Args.OutUnitsDeg logical  = true;
    end
    RAD = 180./pi;
    
    DistAU = convert.length(Args.DistUnits,'au',Dist);
        
    % rectangular ecliptic coordinates of Earth with equinox of J2000
    [E_H] = celestial.SolarSys.calc_vsop87(JD, 'Earth', 'a', 'd');

    % add anti-sun direction at some distance
    SunEarthUnitVector = E_H./sqrt(sum(E_H.^2,1));
    % The Barycentric position of the Earth Shadow at some dist
    E_Sh = E_H + SunEarthUnitVector.*DistAU;
    
     
    Gau = celestial.coo.topocentricVector(JD, Args.GeoPos, 'OutUnits','au',...
                                                             'RefEllipsoid',Args.RefEllipsoid,...
                                                             'Convert2ecliptic',true,...
                                                             'Equinox','J2000');

    % Observer Barycentric position:
    E_H = E_H + Gau;    
    
    % Shadow topocentric position
    T_Sh = E_Sh - E_H;
    
    % Convert to RA/Dec
    
    RA  = atan2(T_Sh(2,:), T_Sh(1,:));
    Dec = atan(T_Sh(3,:)./sqrt( T_Sh(1,:).^2 + T_Sh(2,:).^2  ));
    RA = mod(RA, 2.*pi);

    % no minus sign to E_H because this is already the antiSun direction
    RAantiSun  = atan2(E_H(2,:), E_H(1,:));
    DecantiSun = atan(E_H(3,:)./sqrt( E_H(1,:).^2 + E_H(2,:).^2  ));
    RAantiSun  = mod(RAantiSun, 2.*pi);
    
    if Args.OutUnitsDeg
        RA  = RA.*RAD;
        Dec = Dec.*RAD;
        RAantiSun  = RAantiSun.*RAD;
        DecantiSun = DecantiSun.*RAD;
    end

end