function Data=refellipsoid(RefEllips)
% Return data for a given reference ellipsoid of Earth.
% Package: celestial.Earth
% Description: Return data for a given reference ellipsoid of Earth.
% Input  : - Reference ellipsoid:
%            'Merit1983' - Merit 1983           a=6378137  1/f=298.257
%            'GRS80'     - GRS 80 (IUGG 1980)     6378137      298.257222
%            'GRS67'     - GRS 67 (IUGG 1967)     6378160      298.247167
%            'IAU1976'   - IAU 1976               6378140      298.257
%            'IAU1964'   - IAU 1964               6378160      298.25
%            'WGS84'     - WGS 1984 (default)     6378137      298.257223563
%            'IERS1989'  - IERS 1989              6378136      298.257
% Output : - Data matrix:
%            [Equatorial radius (meters),
%             Flattening factor].
% Reference : Astronomical Almnach
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jun 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Data=celestial.Earth.refellipsoid('WGS84')
% Reliable: 2
%--------------------------------------------------------------------------

switch RefEllips
case {'WGS84'}
    % WGS 1984
    A = 6378137;
    F = 1./298.257223563;
 case {'IERS1989'}
    % IERS 1989
    A = 6378136;
    F = 1./298.257;    
 case {'Merit1983'}
    %'Merit1983' - Merit 1983           a=6378137  1/f=298.257
    A = 6378137;
    F = 1./298.257;
 case {'GRS80'}
    %'GRS80'     - GRS 80 (IUGG 1980)     6378137      298.257222
    A = 6378137;
    F = 1./298.257222;
 case {'GRS67'}
    %'GRS67'     - GRS 67 (IUGG 1967)     6378160      298.247167
    A = 6378160;
    F = 1./298.247167;
 case {'IAU1976'}
    %'IAU1976'   - IAU 1976               6378140      298.257
    A = 6378140;
    F = 1./298.257;
 case {'IAU1964'}
    %'IAU1964'   - IAU 1964               6378160      298.25
    A = 6378160;
    F = 1./298.25;
    
 otherwise
    error('Illegal type for reference ellipsoid');
end
 
Data = [A, F]; 
