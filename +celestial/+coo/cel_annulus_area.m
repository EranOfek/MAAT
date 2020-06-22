function Area=cel_annulus_area(RadIn,RadOut,Units)
% Area within a celestial annulus
% Package: celestial.coo
% Description: Calculate the area within a celestial annulus defined by an
%              inner radius and outer radius of two concentric small circles.
% Input  : - List of Inner radii [radians].
%          - List of Outer radii [radians].
%          - Output Units:
%            'sr'   - sterradians (default).
%            'deg'  - deg^2
%            'am'   - arcmin^2
%            'as'   - arcsec^2
% Output : - List of areas for each pole and radii in output units.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jan 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Area=celestial.coo.cel_annulus_area(0,pi,'deg'); % area of the celestial sphere
% Reliable: 2
%------------------------------------------------------------------------------
RAD        = 180./pi;
DEG_ARCMIN = 60;
DEG_ARCSEC = 3600;

SR_2_Deg2  = RAD.^2;  % conversion factor: sterradian to deg^2

if (nargin==2)
   Units = 'sr';
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end

% calculate sphere cap area [sr]
CapOut = 2.*pi.*(1-cos(RadOut));
CapIn  = 2.*pi.*(1-cos(RadIn));

Area   = CapOut - CapIn;    % [sr]

% convert to output units
switch Units
 case 'sr'
    % do nothing
 case 'deg'
    Area = Area.*SR_2_Deg2;
 case 'am'
    Area = Area.*SR_2_Deg2.*DEG_ARCMIN.^2;
 case 'as'
    Area = Area.*SR_2_Deg2.*DEG_ARCSEC.^2;
 otherwise
    error('Unknown Units Option');
end


