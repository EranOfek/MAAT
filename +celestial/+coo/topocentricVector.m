function [G, Gdot] = topocentricVector(JD_UT1, GeoPos, Args)
    % Calculate the topocentric vector of an observer.
    % Input  : - JD in UT1 time system.
    %          - Geodetic position. If [], then assume geocentric position
    %            and return zeros. Otherwise should be [Long, Lat, Height]
    %            in [rad, rad, m].
    %          * ...,key,val,...
    %            'RefEllipsoid' - Reference ellipsoid. Default is 'WGS84'.
    %            'Convert2ecliptic' - A logical indicating if to convert
    %                   the results to eclitic coordinates.
    %                   Otherwise, the output is in equatorial coordinates.
    %            'Equinox' - Equinox of output: 'date' | 'J2000'.
    %                   Default is 'date'.
    %            'OutUnits' - Output units. Default is 'm'.
    %            'Xp' - The angle of the celestial epheerius pole of the
    %                   Earth with respect to the terrestial pole [rad].
    %                   Along longitude 0. Default is [].
    %            'Yp' - Like Xp but for long of 270 (East). Default is [].
    % Output : - The position vector of the topocentric observer relative
    %            to the Earth center.
    %          - The radius vector time derivative [rad/s].
    % Author : Eran Ofek (Sep 2021)
    % Example: G = celestial.coo.topocentricVector(celestial.time.julday([21 3 2000]), [35 32 0]./RAD)
    
    arguments
        JD_UT1
        GeoPos                              = [];  % [] - geocentric ; [rad, rad, m] - topocentric
        Args.RefEllipsoid                   = 'WGS84';
        Args.Convert2ecliptic(1,1) logical  = false; 
        Args.Equinox                        = 'date';  % 'date' | 'J2000'
        Args.OutUnits                       = 'm';
        Args.Xp                             = [];
        Args.Yp                             = [];
    end
   
    if isempty(GeoPos)
        G  = zeros(3,1);
        Gt = zeros(3,1);
    else
    
        N    = numel(JD_UT1);
        G    = zeros(3,N);
        Gdot = zeros(3,N);
        for I=1:1:N
            
            [~,GeocCart] = celestial.Earth.geod2geoc(GeoPos, Args.RefEllipsoid);

            LAST = celestial.time.lst(JD_UT1(I), 0, 'a');     % calculate app. sidereal time at Greenwich
            LAST = LAST.*2.*pi;             % convert to radians

            RotMat = tools.math.geometry.rotm(LAST,3);   % and not -LAST !!!
            if ~isempty(Args.Xp) && ~isempty(Args.Yp)
                RotXY = tools.math.geometry.rotm(Args.Yp,1)*tools.math.geometry.rotm(Args.Xp,2);
            else
                RotXY = diag(ones(1,3));
            end

            G(:,I) = RotMat*RotXY*GeocCart.';
        end
        
        if Args.Convert2ecliptic
            RotMatEq2Ec = celestial.coo.rotm_coo('e');
            G= RotMatEq2Ec * G;
        end

        switch lower(Args.Equinox)
            case 'date'
                % do nothing
            case 'j2000'
                RotP = celestial.coo.rotm_coo('pd', 2451545.5);
                G    = RotP * G;
            otherwise
                error('Unknown Equinox option');
        end

        if nargout>1
            W = 7.2921151467e-5;  % [rad/s]
            RotRot = [-sin(LAST), -cos(LAST), 0; cos(LAST), -sin(LAST), 0; 0, 0, 0];
            Gdot = W.*RotRot.*RotXY*GeocCart;
            if Args.Convert2ecliptic
                G = RotMatEq2Ec * G;
            end

            switch lower(Args.Equinox)
                case 'date'
                    % do nothing
                case 'j2000'
                    RotP = celestial.coo.rotm_coo('pd');
                    Gdot    = RotP * Gdot;
                otherwise
                    error('Unknown Equinox option');
            end
        end

        G = convert.length('m',Args.OutUnits, G);
    end
end

