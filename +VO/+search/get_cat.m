function Cat=get_cat(CatName,RA,Dec,Radius,varargin)
% Search selected astronomical catalogs
% Package: VO.search
% Description: Search selected astronomical catalogs
% Input  : - Catalog name: 'sdss' | 'apass' | ...
%          - J2000.0 R.A. [default radians or sexagesimal].
%          - J2000.0 Dec. [default radians or sexagesimal].
%          - Radius [default in arcsec].
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'CooUnits' - input coordinate units 'rad'|'deg'.
%                         Default is 'rad'.
%            'RadUnits' - 'arcsec'|'arcmin'|'deg'|'rad'.
%                         Default is 'arcsec'.
%            'Shape'    - 'circ'. Default is 'circ'.
% Output : - An AstCat object with the catalog
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=VO.search.get_cat('sdss',pi,0.1,3600);
% Reliable: 2
%--------------------------------------------------------------------------
RAD = 180./pi;


DefV.CooUnits             = 'rad';
DefV.RadUnits             = 'arcsec';
DefV.Shape                = 'circ';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (ischar(RA))
    RA = celestial.coo.convertdms(RA,'SH','r');
    InPar.CooUnits = 'rad';
end
if (ischar(Dec))
    Dec = celestial.coo.convertdms(Dec,'SD','R');
    InPar.CooUnits = 'rad';
end

RA     = convert.angular(InPar.CooUnits,'rad',RA);
Dec    = convert.angular(InPar.CooUnits,'rad',Dec);
Radius = convert.angular(InPar.RadUnits,'arcsec',Radius);

switch lower(CatName)
    case 'sdss'
        Cat = VO.search.get_sdss(RA,Dec,Radius);
        
    case 'apass'
        Cat = VO.search.get_apass(RA,Dec,Radius./(3600.*RAD));
        
    otherwise
        error('Unknown CatName option');
end


        
        