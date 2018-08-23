function [PA]=parallactic_angle(varargin)
%------------------------------------------------------------------------------
% parallactic_angle function                                             ephem
% Description: Calculate the parallactic angle of an object.
%              The parallactic is defined as the angle between the local
%              zenith, the object and the celestial north pole measured
%              westwerd (e.g., negative before, and positive after the
%              passage through the southern meridian).
% Input  : * Set of three input argument: [RA, Dec], LST, Lat
%            or alternatively four input argument: RA, Dec, LST, Lat.
%            Where RA, Dec and Lat are in radians, and LST in 
%            fraction of days. Lat is observer the geodetic latitude.
%            LST can be either a scalar, matrix of the same size as
%            RA and Dec, or a vector which have a common dimension
%            as RA and Dec.
% Output : - Parallactic angle
%            If object in the zenith then NaN.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                   October 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 1
%------------------------------------------------------------------------------

if (length(varargin)==3),
   Coo = varargin{1};
   LST = varargin{2};
   Lat = varargin{3};

   RA  = Coo(:,1);
   Dec = Coo(:,2);
else
   RA   = varargin{1};
   Dec  = varargin{2};
   LST  = varargin{3};
   Lat  = varargin{4};
end
%HA = 2.*pi.*LST - RA;
HA = bsxfun(@minus,2.*pi.*LST,RA);

TanQ = sin(HA)./(tan(Lat).*cos(Dec) - sin(Dec).*cos(HA));
PA   = atan(TanQ);
