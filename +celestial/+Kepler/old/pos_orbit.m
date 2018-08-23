function [EqCoo0,AngSpeed,Mag]=pos_orbit(Time,OrbEl,EarthPosFun,MagFun,varargin);
%--------------------------------------------------------------------------
% pos_orbit function                                                 ephem
% Description: Given orbital elements, calculate the Geocentric position
%              of the celestial body. The program can calculate the
%              positions at a single time for many bodies or the
%              positions at multiple times for a single body.
% Input  : - Time in TT time scale on which to calculate bodies position.
%            This can be in JD, or [D M Y H M S] or [D M Y Frac],...
%            If scalar than, calculte positions for multiple bodies.
%            If vector/matrix than assumes a single body.
%          - Orbital elements [J2000.0; radians/au]
%            matrix of the form: [T, q, e, LongPeri, LongAsc, Inc]
%            or a structure containing the appropriate fields
%            from the following list:
%              .PeriJD
%              .Q
%              .E
%              .LonPeri
%              .LonAsc
%              .Inc
%          - Earth position function:
%            {'vsop87' | 'low'}, default is 'low'.
%          - Magnitude model function:
%            Examples: 
%          * Arbitrary number of input arguments to be passed to the
%            magnitude model function.
%            Examples:
% Output : -
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Jan 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%--------------------------------------------------------------------------
RAD = 180./pi;

if (size(Time,2)>1),
   JD = julday(Time).';
else
   JD = Time;
end


% Read orbital elements
if (isstruct(OrbEl)==0),
   El.PeriJD = OrbEl(:,1);
   El.Q      = OrbEl(:,2);
   El.E      = OrbEl(:,3);
   El.LonPeri= OrbEl(:,4);
   El.LonAsc = OrbEl(:,5);
   El.Inc    = OrbEl(:,6);   
else
   El = OrbEl;
end

if (length(JD)==1),
   %-------------------------------------
   %--- single time; multiple objects ---
   %-------------------------------------

   [Nu,R,E,Vel] = kepler_elliptic(JD-El.PeriJD,El.Q,El.E);
   XYZ = trueanom2pos(R,Nu,El.LonAsc,El.LonPeri,El.Inc);
   [Ex,Ey,Ez]=ple_planet(JD,'Earth','XYZ');   % equinox of date!
   X0 = XYZ(:,1) - Ex;
   Y0 = XYZ(:,2) - Ey;
   Z0 = XYZ(:,3) - Ez;
   Long0   = atan2(Y0,X0);
   Lat0    = atan(Z0./sqrt(X0.^2 + Y0.^2));
   EqCoo0  = coco([Long0, Lat0],'e','j2000.0');

   R = sqrt(sum(XYZ.^2,2));
   Delta = sqrt(X0.^2 + Y0.^2 + Z0.^2);

   Nobj = length(R);
   V1 = zeros(Nobj,1);
   V2 = zeros(Nobj,1);
   for J=1:1:Nobj
      W1 = XYZ(J,:);
      V1 = W1./norm(W1);
      W2 = [X0(J), Y0(J), Z0(J)];
      V2 = W2./norm(W2);

      PhaseAng = acos(dot(V1,V2));
   end

   DelT = 0.01;
   [Nu,R,E,Vel] = kepler_elliptic(JD+DelT-El.PeriJD,El.Q,El.E);
   XYZ = trueanom2pos(R,Nu,El.LonAsc,El.LonPeri,El.Inc);
   [Ex,Ey,Ez]=ple_planet(JD+DelT,'Earth','XYZ');  % equinox of date!
   X1 = XYZ(:,1) - Ex;
   Y1 = XYZ(:,2) - Ey;
   Z1 = XYZ(:,3) - Ez;
   Long1   = atan2(Y1,X1);
   Lat1    = atan(Z1./sqrt(X1.^2 + Y1.^2));
   EqCoo1  = coco([Long1, Lat1],'e','j2000.0');

   D = sphere_dist(EqCoo0(:,1),EqCoo0(:,2), EqCoo1(:,1),EqCoo1(:,2));
   AngSpeed = D./DelT;    % rad/day
   Mag = asteroid_magnitude(R,Delta,PhaseAng,varargin{1},varargin{2});

else
   %------------------------------------
   %--- multiple time; single object ---
   %------------------------------------



end

