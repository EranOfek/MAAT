function [BestCen,BestRad,BestRMS,BestDif,BestRDi,Info]=fit_circle(Data,Geom,MinPar,Nstep);
%------------------------------------------------------------------------------
% fit_circle function                                                   FitFun
% Description: Fit points, on a plane or a sphere, to a circle.
%              Calculate the best fit radius and the center of the circle.
% Input  : - Matrix of data [RA, Dec, Property].
%            Where RA can be a column vector [radians],
%            or matrix of [H M S].
%            Dec can be a column vector [radians],
%            or matrix of [Sign D M S].
%            An optional column of a property of each point.
%          - Geometry:
%            'sphere' - celestial sphere.
%            'plane'  - plane geometry.
%          - Minimizing method:
%            'rms'  - by radii rms
%            'diff' - by max(R)-min(R)
%            'rd'   - by (max(R)-min(R))/R
%          - Number of scanning steps, default is 200.
% Output : - Best center [X,Y] or [RA,Dec] in radians.
%          - Best radius in radian.
%          - Best RMS
%          - Best max(R)-min(R)
%          - Best (max(R)-min(R))/R
%          - Information on each point, one line per point, sorted by PA:
%            [Distance from best center [radians],
%             PA relative to best center [radians],
%             Property of point (if not given then = NaN),
%             Clockwise angular distance from next point on the circle [rad]
%             Index of respective line in original data matrix]
% Tested : Matlab 5.3
%     By : Eran O. Ofek                   January 2004
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 1
%------------------------------------------------------------------------------

RAD = 180./pi;

if (nargin==3),
   Nstep = 200;
elseif (nargin==4),
   % do nothing
else
   error('Illigal number of inpur arguments');
end

LD = length(Data(1,:));

if (LD==2 | LD==3),
   RA  = Data(:,1);
   Dec = Data(:,2);
   if (LD==3),
      Prop = Data(:,3);
   else
      Prop = NaN.*RA;
   end
elseif (LD==7 | LD==8),
   RA  = convertdms(Data(:,1:3),'H','r');
   Dec = convertdms(Data(:,4:7),'D','R');
   if (LD==8),
      Prop = Data(:,8);
   else
      Prop = NaN.*RA;
   end
else
   error('Illigal number of columns in Data');
end

Np      = length(RA);

MeanRA  = mean(RA);
MeanDec = mean(Dec);
switch Geom
 case 'sphere'
    DistIJ  = celestial.coo.sphere_dist(RA(1),Dec(1),RA(2:end),Dec(2:end));
 case 'plane'
    DistIJ  = Util.Geom.plane_dist(RA(1),Dec(1),RA(2:end),Dec(2:end));
 otherwise
    error('Unknown geometry option');
end


MaxDist = max(DistIJ);

Step     = 2.*MaxDist./Nstep;
X        = [MeanRA-MaxDist:Step:MeanRA+MaxDist];
Y        = [MeanDec-MaxDist:Step:MeanDec+MaxDist];
Nx       = length(X);
Ny       = length(Y);

Xmat     = ones(Ny,1)*X;
Ymat     = Y.'*ones(1,Nx);

Dist     = zeros(Ny,Nx,Np);
switch Geom
 case 'sphere'
    for I=1:1:Np
       DX            = Xmat - RA(I);
       Dist(:,:,I)   = acos(sin(Ymat).*sin(Dec(I)) + cos(Ymat).*cos(Dec(I)).*cos(DX));
    end
 case 'plane'
    for I=1:1:Np
       Dist(:,:,I)   = sqrt((Xmat-RA(I)).^2 + (Ymat-Dec(I)).^2);
    end
 otherwise
    error('Unknown geometry option');
end

switch MinPar
 case 'rms'
    MeanDist     = mean(Dist,3);
    RMS          = sqrt(mean((Dist-repmat(MeanDist,[1 1 Np])).^2));
    [Min,MinInd] = min2d(RMS);
 case 'diff' 
    Diff         = max(Dist,[],3) - min(Dist,[],3);
    [Min,MinInd] = min2d(Diff);
 case 'rd'
    MeanDist     = mean(Dist,3);
    RD           = (max(Dist,[],3) - min(Dist,[],3))./MeanDist;   
    [Min,MinInd] = min2d(RD);
 otherwise
    error('Unknown MinPar option');
end

I = MinInd(1);
J = MinInd(2);
BestCen   = [X(J),Y(I)];
switch Geom
 case 'sphere'
    [Dist,PA] = celestial.coo.sphere_dist(RA,Dec,BestCen(1),BestCen(2)); 
 case 'plane'
    [Dist,PA] = Util.Geom.plane_dist(RA,Dec,BestCen(1),BestCen(2)); 
 otherwise
    error('Unknown geometry option');
end

BestRad = mean(Dist);
BestRMS = sqrt(mean((Dist - BestRad).^2));
BestDif = max(Dist) - min(Dist);
BestRDi = (max(Dist) - min(Dist))./BestRad;
BestDPA = [Dist, PA, Prop];

%BestRMS.*RAD.*3600
%BestCen.*RAD
%BestRad.*RAD.*3600
%BestDif.*RAD.*3600
%BestRDi


%--- Distance to next (counter clockwise) point on circle (deg) ---

[Info,IndInfo] = sortrows(BestDPA,2);
Ni             = length(Info(:,1));
DiffPA         = [Info(1,2)+2.*pi-Info(Ni,2) ; diff(Info(:,2))];
Info           = [Info, DiffPA, IndInfo];

