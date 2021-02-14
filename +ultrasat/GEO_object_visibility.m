function GeoObjVis=GEO_object_visibility(JD,Coo,varargin)
% SHORT DESCRIPTION HERE
% Package: ultrasat
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings
% Output : - 
% License: GNU general public license version 3
%     By : Yossi Shvartzvald              Jan 2021
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: N1  = [220./RAD, 66./RAD];
%          S1  = [ 42./RAD,-66./RAD];
%          Coo = [N1;S1];
%          JD  = celestial.time.julday([1 1 2025 0]) + (0:0.1:365)';
%          GeoObjVis=ultrasat.GEO_object_visibility(JD,Coo);
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;
D_MOON = 384400;
D_GEO  = 36000;
EarthSD = 86164.091./86400;

DefV.Name                 = {};     % override all the parameters provided in the list
DefV.JD_vernal_equinox = celestial.time.julday([21 3 2023 0]);
DefV.GEO_slot = -4/RAD; %4W


if (numel(varargin)==1)
    % assume input is a structure (like DefV)
    InPar = varargin{1};
else
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);
end

if (~isempty(InPar.Name))
    Pars = DefPar.(InPar.Name);
    for I=1:2:(numel(Pars)-1)
        InPar.(Pars{I}) = Pars{I+1};
    end
end

[SunRA,SunDec]     = celestial.SolarSys.suncoo(JD,'j');

EarthRA   = 2.*pi.*mod(JD-InPar.JD_vernal_equinox,EarthSD)./EarthSD+InPar.GEO_slot;
EarthDec  = zeros(numel(JD),1);        %assuming Equatorial orbit

[MoonRA,MoonDec]   = celestial.SolarSys.mooncool(JD,NaN,'l');
[xM,yM,zM] = sph2cart(MoonRA,MoonDec,D_MOON);
[xE,yE,zE] = sph2cart(EarthRA,EarthDec,D_GEO);
[MoonRA,MoonDec,~]= cart2sph(xM+xE,yM+yE,zM+zE);

GeoObjVis.SunAngDist = celestial.coo.sphere_dist_fast(Coo(:,1)',Coo(:,2)',SunRA,SunDec);
GeoObjVis.EarthAngDist = celestial.coo.sphere_dist_fast(Coo(:,1)',Coo(:,2)',EarthRA,EarthDec);
GeoObjVis.MoonAngDist = celestial.coo.sphere_dist_fast(Coo(:,1)',Coo(:,2)',MoonRA,MoonDec);

GeoObjVis.JD = JD;
GeoObjVis.Coo = Coo;




