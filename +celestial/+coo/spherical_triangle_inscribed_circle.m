function r=spherical_triangle_inscribed_circle(Method,varargin)
% Calculate the radius of the inscribed circle of a spherical triangle
% Package: celestial.coo
% Input  : - Method by which to calculate the radius. The following methods
%            are supported:
%            'sides' - calculate the radius from the sides a,b,c.
%            'angles'- calculate the radius from the angles A,B,C.
%            'vertex'- calculate the radius from the triangle vertces.
%          * 3 or 6 arguments describing the spherical triangles.
%            For 'sides', the arguments are arrays of a,b,c.
%            For 'angles', the arguments are arrays of A,B,C.
%            For 'vertex', the arguments are arrays of
%            Long1,Lat1,Long2,Lat2,Long3,Lat3.
%            All arguments are in radians.
% Output : - An array of the radii of the inscribed circles.
% Reference: Todhunter & Leathem 
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: r=celestial.coo.spherical_triangle_inscribed_circle('sides',1,1,1)
% Reliable: 2


import celestial.coo.*

switch lower(Method)
    case 'sides'
        % a,b,c
        [a,b,c] = deal(varargin{:});
        
        s = 0.5.*(a + b + c);
        r = atan(sqrt( sin(s-a).*sin(s-b).*sin(s-c)./sin(s) ));
    case 'angles'
        % A,B,C
        [A,B,C] = deal(varargin{:});
        
        S = 0.5.*(A + B + C);
        r = atan( sqrt(-cos(S).*cos(S-A).*cos(S-B).*cos(S-C) )./(2.*cos(0.5.*A).*cos(0.5.*B).*cos(0.5.*C) ) );
    case 'vertex'
        [Long1,Lat1,Long2,Lat2,Long3,Lat3] = deal(varargin{:});        a = sphere_dist(Long2,Lat2,Long3,Lat3);
        b = sphere_dist(Long3,Lat3,Long1,Lat1);
        c = sphere_dist(Long1,Lat1,Long2,Lat2);
        
        s = 0.5.*(a + b + c);
        r = atan(sqrt( sin(s-a).*sin(s-b).*sin(s-c)./sin(s) ));
    otherwise
        error('Unknown Method option')
end