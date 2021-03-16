function U2=aberration(U1,V,Units)
% Apply aberration of light to source position
% Package: celestial.coo
% Description: Calculate the position of a star corrected for aberration
%              of light.
%              Rigoursly, these is applied after accounting for the light
%              travel time effect and light deflection and before the
%              precession and nutation are being applied.
% Input  : - Matrix containing three columns (cosine directions).
%            Each row represent the object position in a given epoch.
%            If two columns are given the assume J2000.0 [RA, Dec] in radians.
%          - The observer velocity in the barycentric reference frame,
%            in the same format as the previous parameter (Equatorial J2000.0).
%            If only one column is specified, then this is assumed
%            to be the JD and Geocentric observer. In this case the progarm
%            will calculate the Geocentric velocity vectors for each epoch
%            using the Ron & Vondrak (1986) expressions.
%          - Input/output vector units - options are:
%            'cgs'  - for [cm].
%            'SI'   - for [m].
%            'agd'  - for [AU] - default.
% Output : - The position of the object corrected for abberation of light
%            and given in the observer reference frame in the same
%            format as the input parameters.
% Reference: Explanatory supplement to the Astronomical Almanac, p. 150.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Oct 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: U2=celestial.coo.aberration([celestial.coo.convertdms([2 44 12.9747],'H','r'), celestial.coo.convertdms([1 49 13 39.896],'D','R')],celestial.time.julday([13.19 11 2028]));
% Reliable: 2
%------------------------------------------------------------------------------

Def.Units = 'agd';
if (nargin==2)
   Units = Def.Units;
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end

switch lower(Units)
 case 'agd'
    C = constant.c.*convert.units('cm/s','au/day');
 otherwise
    C = constant.c(Units);
end



if (size(U1,2)==3)
   % already in cosine direction
   OutType = 'cosd';
elseif (size(U1,2)==2)
   % convert RA,Dec to cosine directions:
   OutType = 'rad';
   U1 = celestial.coo.cosined(U1);
else
   error('Position matrix has illegal number of columns');
end   

if (size(V,2)==1)
   % calculate the Geocentric velocity using the Ron & Vondrak (1986)
   % expressions
   JD = V;
   V = celestial.SolarSys.earth_vel_ron_vondrak(JD,Units);

elseif (size(V,2)==3)
   % assume each row represent the Geocentric velocity
else
   error('Velocity matrix has illegal number of columns');
end


% new section
Nv = size(V,1);
Nu = size(U1,1);
V       = V./C;
if (Nv==Nu)
   % do nothing
else
   if (Nu>Nv)
      V = repmat(V,Nu,1);
   else
      U1 = repmat(U,Nv,1);
   end
end

AbsU1   = sqrt(sum(U1.^2,2));
NormU1  = bsxfun(@rdivide,U1,AbsU1);
AbsV    = sqrt(sum(V.^2,2));
InvBeta = sqrt(1-AbsV.^2);

DotUV = dot(NormU1,V,2);

U2 = bsxfun(@rdivide,(bsxfun(@times,InvBeta,NormU1) + V + bsxfun(@rdivide,bsxfun(@times,DotUV,V), 1 + InvBeta)),(1 + DotUV));




% old section
if (1==0)

AbsU1   = sqrt(sum(U1.^2,2));
NormU1  = bsxfun(@rdivide,U1,AbsU1);
V       = V./C;
AbsV    = sqrt(sum(V.^2,2));
InvBeta = sqrt(1-AbsV.^2);

F1      = dot(NormU1,V,2);
F2      = 1 + F1./(1 + InvBeta);

U2      = (bsxfun(@times,U1,InvBeta) + ...
           bsxfun(@rdivide,bsxfun(@times,F2.*AbsU1,V),(1 + F1) ));

end


switch lower(OutType)
 case 'cosd'
    % do nothing
 case 'rad'
    % convert cosine directions to [RA, Dec] in radians
    U2 = celestial.coo.cosined(U2);

    Ineg = find(U2(:,1)<0);
    U2(Ineg,1) = 2.*pi + U2(Ineg,1);
 otherwise
    error('Unknown OutType option');
end
