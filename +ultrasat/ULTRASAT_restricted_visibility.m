function ULTRASAT_vis=ULTRASAT_restricted_visibility(JD,Coo,varargin)
% SHORT DESCRIPTION HERE
% Package: ultrasat
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Yossi Shvartzvald                    Jan 2021
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: N1  = [220./RAD, 66./RAD];
%          S1  = [ 42./RAD,-66./RAD];
%          Coo = [N1;S1];
%          JD  = celestial.time.julday([1 1 2025 0]) + (0:0.1:365)';
%          ULTRASAT_vis=ultrasat.ULTRASAT_restricted_visibility(JD,Coo);
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.Name                 = {};     % override all the parameters provided in the list

DefV.MinSunDist = 70./RAD;
DefV.MinMoonDist = 49./RAD;
DefV.MinEarthDist = 56./RAD;

DefV.Power_MaxSunDist = 135./RAD; %90+45; 45 for thr solar panels
DefV.Comm_MinEarthDist = 40./RAD;
DefV.Comm_MaxEarthDist = 140./RAD;

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

ULTRASAT_vis = ultrasat.GEO_object_visibility(JD,Coo);

ULTRASAT_vis.SunLimits   = ULTRASAT_vis.SunAngDist   > InPar.MinSunDist;
ULTRASAT_vis.EarthLimits = ULTRASAT_vis.EarthAngDist > InPar.MinEarthDist;
ULTRASAT_vis.MoonLimits  =ULTRASAT_vis.MoonAngDist  > InPar.MinMoonDist;

ULTRASAT_vis.PowerLimits = ULTRASAT_vis.SunAngDist   < InPar.Power_MaxSunDist;
ULTRASAT_vis.CommLimits  = ULTRASAT_vis.EarthAngDist > InPar.Comm_MinEarthDist & ...
                           ULTRASAT_vis.EarthAngDist < InPar.Comm_MaxEarthDist;