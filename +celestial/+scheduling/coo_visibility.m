function [Res,LimitAM]=coo_visibility(JD,RA,Dec,varargin)
% Calculate the visibility of celestial coordinates
% Package: +celestial.scheduling
% Description: For a given night, and a list of targets, calculate the
%              matrix of target visibility in 5-min (default) steps.
%              The visibility criteria, includes SunAltLimit, MaxAM,
%              Moon distance as a function of illumination, and minimum
%              visibility time. By default it will also calculate the AM
%              limit above the target is found for
%              MinVisibilityTime*FactorVisibilityAM.
% Input  : - JD. will choose the nearest night (unless AllNightVisibility
%            is false).
%          - Vector of targets J2000.0 RA [rad].
%          - Vector of targets J2000.0 Dec [rad].
%          * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'AllNightVisibility' - If true, then calculate the all night
%                   visibility of the targets 9e.g., every TimeStep).
%                   If false, then calculate the visibility only at the
%                   current requested time.
%                   Default is true.
%            'Lon' - Observatory East longitude [rad]. Default is 35./RAD.
%            'Lat' - Observatory North Latitude [rad]. Default is 30./RAD.
%            'TimeStep' - Time step [day] for time visibility table.
%                   Default is 5./1440.
%            'SunAltLimit' - Sun altitude limit [rad]. Default is -12./RAD.
%            'MaxAM' - Maximum AM limit. Default is 2.
%            'AzAltLimit' - Table of allowed altitude vs. Az [deg].
%                   Default is [0 15;90 15; 180 15; 270 15; 360 15]
%                   Must start at Az 0 and end at Az 360.
%            'MinMoonDistIllum' - Minimum angular distance from the moon
%                   [deg] as a afunction of Moon illumination.
%                   Default is
%                   [0 0; 0.1 1; 0.2 1; 0.3 1; 0.4 2; 0.5 3; 0.6 5;0.7 10;0.8 15; 0.9 30; 1.0 30]
%            'MinVisibilityTime' - Minimum required target time visibility
%                   during the night [day]. Default is 1./24.
%            'FactorVisibilityAM' - Factor by which to multiply the
%                   visibility time after adjusting it to new AM.
%            'CleanTime' - Remove times in which Sun is above alt limit.
%                   Default is true.
%            'InterpMethod' - Tables interpolation. Default is 'linear'
% Output : - A structure with the following fields:
%            .VecJD - A vector of JD for each line in the visibility
%                   tables.
%            .Flag - A matrix of logicals (Time,Target) indicating if the
%                   target is visible in each time step.
%            .RA - Vector of J2000.0 RA [rad].
%            .Dec - Vector of J2000.0 Dec [rad].
%            .TimeVisible - Vector indicating the total time during the
%                   night in which the target is visible [day].
%            .NtargetVisible - Vector of te number of targets visible on
%                   each time.
%            .AM - Matrix of AM for each target and time.
%            .LST - LST [rad]
%            .MoonDist - Matrix of moon distances [rad].
%            .FlagLim - The same as Flag, but after applying the new AM
%                   limit per target.
%            .TimeVisibleLim - The same as TimeVisible, but after applying the new AM
%                   limit per target.
%            .NtargetVisibleLim - The same as NtargetVisible, but after applying the new AM
%                   limit per target.
%          - A vector of AM limits per target as calculated by requireing
%            the MinVisibilityTime for the target
%      By: Eran Ofek                      Oct 2020
% Example: [TileList,TileArea]=celestial.coo.tile_the_sky(56,42);
% [Res,LimitAM]=celestial.scheduling.coo_visibility(2451556,TileList(:,1),TileList(:,2))
% [Res]=celestial.scheduling.coo_visibility(2451556,TileList(:,1),TileList(:,2),'MaxAM',LimitAM)


RAD = 180./pi;

InPar = inputParser;
addOptional(InPar,'AllNightVisibility',true);
addOptional(InPar,'Lon',35./RAD);
addOptional(InPar,'Lat',30./RAD);
addOptional(InPar,'TimeStep',5./1440);  % day
addOptional(InPar,'SunAltLimit',-12./RAD);  
addOptional(InPar,'MaxAM',2);  
addOptional(InPar,'AzAltLimit',[0 15;90 15; 180 15; 270 15; 360 15]);  % [Az Alt] deg  
addOptional(InPar,'MinMoonDistIllum',[0 0; 0.1 1; 0.2 1; 0.3 1; 0.4 2; 0.5 3; 0.6 5;0.7 10;0.8 15; 0.9 30; 1.0 30]);  % [illum, MinDist]  
addOptional(InPar,'MinVisibilityTime',1./24);  % [day] 
addOptional(InPar,'FactorVisibilityAM',1.2);  % [day] 
addOptional(InPar,'CleanTime',true); % remove times in which the Sun is below limit
addOptional(InPar,'InterpMethod','linear'); 

parse(InPar,varargin{:});
InPar = InPar.Results;

if InPar.AzAltLimit(1,1)~=0 || InPar.AzAltLimit(end,1)~=360
    error('AzAltLimit must contain info for Az 0 and 360 deg');
end

MinNptVisibility = InPar.MinVisibilityTime./InPar.TimeStep;

% objects are in columns
RA  = RA(:).';
Dec = Dec(:).';
InPar.MaxAM = InPar.MaxAM(:).';

if InPar.AllNightVisibility
    % time coordinate is in rows
    FloorJD = floor(JD);
    MidDayJD = FloorJD - InPar.Lon./(2.*pi);
    VecJD = (MidDayJD:InPar.TimeStep:(MidDayJD+1)).';
else
    VecJD = JD;
end

[SunRA,SunDec]     = celestial.SolarSys.suncoo(VecJD,'j');

LST                = celestial.time.lst(VecJD,InPar.Lon,'m').*2.*pi;  % [rad]

[SunAz,SunAlt]     = celestial.coo.hadec2azalt(LST-SunRA,SunDec,InPar.Lat);

if InPar.CleanTime
    FlagSun = SunAlt<InPar.SunAltLimit;
    
    VecJD  = VecJD(FlagSun);
    SunRA  = SunRA(FlagSun);
    SunDec = SunDec(FlagSun);
    SunAz  = SunAz(FlagSun);
    SunAlt = SunAlt(FlagSun);
    LST    = LST(FlagSun);
end

[MoonRA,MoonDec]   = celestial.SolarSys.mooncool(VecJD,[InPar.Lon, InPar.Lat],'b');
[MoonIllum,MoonPh] = celestial.SolarSys.moon_illum(VecJD);
    
[MoonAz,MoonAlt]   = celestial.coo.hadec2azalt(LST-MoonRA,MoonDec,InPar.Lat);
MoonDist           = celestial.coo.sphere_dist_fast(RA,Dec,MoonRA,MoonDec);

MinMoonDist = interp1(InPar.MinMoonDistIllum(:,1),InPar.MinMoonDistIllum(:,2),abs(MoonIllum),InPar.InterpMethod)./RAD;  % [rad]

% HA is a matrix of [time, obj]
HA                 = LST - RA;
[Az,Alt]           = celestial.coo.hadec2azalt(HA,Dec,InPar.Lat);
AM                 = celestial.coo.hardie(pi./2-Alt);

AltLimit           = interp1(InPar.AzAltLimit(:,1)./RAD,InPar.AzAltLimit(:,2)./RAD,Az,InPar.InterpMethod);


Flag = AM<InPar.MaxAM & Alt>AltLimit & SunAlt<InPar.SunAltLimit & MoonDist>MinMoonDist;
Flag = Flag & sum(Flag)>MinNptVisibility;

Res.VecJD          = VecJD;
Res.Flag           = Flag;
Res.RA             = RA;
Res.Dec            = Dec;
Res.TimeVisible    = sum(Res.Flag).*InPar.TimeStep;  % [day]
Res.NtargetVisible = sum(Flag,2);
Res.AM             = AM;
Res.MoonDist       = MoonDist;
Res.LST            = LST;

if InPar.AllNightVisibility && nargout>1
    AM(AM>InPar.MaxAM) = NaN;
    % check what is the airmass below, which the object reside at least
    % InPar.MinVisibilityTime*InPar.FactorVisibilityAM
    SortedAM = sort(AM,1);
    
    RequiredNptBelowAM = ceil(InPar.MinVisibilityTime*InPar.FactorVisibilityAM./InPar.TimeStep);
    
    LimitAM = SortedAM(RequiredNptBelowAM,:);
    LimitAM(isnan(LimitAM)) = InPar.MaxAM;
    
    
    FlagLim = AM<LimitAM & Alt>AltLimit & SunAlt<InPar.SunAltLimit & MoonDist>MinMoonDist;
    FlagLim = FlagLim & sum(FlagLim)>MinNptVisibility;

    Res.FlagLim           = FlagLim;
    Res.TimeVisibleLim    = sum(Res.FlagLim).*InPar.TimeStep;  % [day]
    Res.NtargetVisibleLim = sum(FlagLim,2);
    
else
    LimitAM = [];
end