function [R,Long,Lat]=spherical_triangle_circum_circle(Method,varargin)
% Calculate the radius of the circum circle of a spherical triangle
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
% Output : - An array of the radii of the circum circles.
% Reference: Todhunter & Leathem 
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: R=celestial.coo.spherical_triangle_circum_circle('sides',1,1,1)
% Reliable: 2


import celestial.coo.*

switch lower(Method)
    case {'s','sides'}
        % a,b,c
        [a,b,c] = deal(varargin{:});
        
        s = 0.5.*(a + b + c);
        R = atan(2.*sin(0.5.*a).*sin(0.5.*b).*sin(0.5.*c)./sqrt(sin(s).*sin(s-a).*sin(s-b).*sin(s-c) ));
    case {'a','angles'}
        % A,B,C
        [A,B,C] = deal(varargin{:});
        
        S = 0.5.*(A + B + C);
        R = atan(sqrt( -cos(S)./(cos(S-A).*cos(S-B).*cos(S-C) ) ) );
    case {'v','vertex'}
        [Long1,Lat1,Long2,Lat2,Long3,Lat3] = deal(varargin{:});       
        [a,PAa] = sphere_dist(Long2,Lat2,Long3,Lat3);
        [b,PAb] = sphere_dist(Long3,Lat3,Long1,Lat1);
        [c,PAc] = sphere_dist(Long1,Lat1,Long2,Lat2);
        
        s = 0.5.*(a + b + c);
        R = atan(2.*sin(0.5.*a).*sin(0.5.*b).*sin(0.5.*c)./sqrt(sin(s).*sin(s-a).*sin(s-b).*sin(s-c) ));
    otherwise
        error('Unknown Method option')
end

if (nargout>1)
    switch lower(Method)
        case {'v','vertex'}
            % calculate the center of the circum circle
            % PA is measure relative to [Long1,Lat1] in the direction
            % between p2 and p3
            %PAc, PAb
            
             PA = acos(tan(0.5.*c).*cot(R))
%             Az = PAc +/- PA;
%             [Lat,Long]=reckon(Lat1,Long1,R,Az,'radians');
        otherwise
            error('2nd and 3rd output are allowed only with the Method=vertex option');
    end
end