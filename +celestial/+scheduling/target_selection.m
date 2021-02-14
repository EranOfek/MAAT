function ResS=target_selection(varargin)
% Select a single target for one telescope 
% 
% Package: +celestial.scheduling
% Example: Res=celestial.scheduling.target_selection

RAD = 180./pi;

InPar = inputParser;
addOptional(InPar,'JD',[]); % for which to do the scheduling
addOptional(InPar,'TargetList',[]); % [RA, Dec, LastObs, MainCounter, NightCounter]


addOptional(InPar,'Visibility',[]); % Structure output from celestial.scheduling.coo_visibility for all targets

addOptional(InPar,'Lon',35./RAD);
addOptional(InPar,'Lat',30./RAD);
addOptional(InPar,'ExpTime',5./1440);


addOptional(InPar,'DecRange',[-30 90]./RAD);


% parameters for visibility calculation
addOptional(InPar,'TimeStep',5./1440);  % day
addOptional(InPar,'SunAltLimit',-12./RAD);  
addOptional(InPar,'MaxAM',2);   % or vector
addOptional(InPar,'AzAltLimit',[0 15;90 15; 180 15; 270 15; 360 15]);  % [Az Alt] deg  
addOptional(InPar,'MinMoonDistIllum',[0 0; 0.1 1; 0.2 1; 0.3 1; 0.4 2; 0.5 3; 0.6 5;0.7 10;0.8 15; 0.9 30; 1.0 30]);  % [illum, MinDist]  
addOptional(InPar,'MinVisibilityTime',5./24);  % [day] 
addOptional(InPar,'FactorVisibilityAM',1.2);  % [day] 

% parameters for weight calculation
addOptional(InPar,'MainCadence',0.8);  % [day]
addOptional(InPar,'NightCadence',30./1440); % [day]
addOptional(InPar,'Nfast',8); % [day]
addOptional(InPar,'MainWFun',@celestial.scheduling.fermiexp); %@(t) 1.0+0.5.*exp(-t./1) ); % weight as a function of time since it is allowed to observe the target
addOptional(InPar,'MainWFunPar', {0.8, 1, 0.03, 1, 0.5} );  %t0, DecayExp, SoftFermi, BaseW, ExtraW
addOptional(InPar,'NightWFun',@celestial.scheduling.fermiexp); %@(t) 1.5+0.5.*exp(-t./1) ); % weight as a function of time since it is allowed to observe the target
addOptional(InPar,'NightWFunPar',{30./1440, 1, 0.003, 1.5, 0.5});  %t0, DecayExp, SoftFermi, BaseW, ExtraW

addOptional(InPar,'InterpMethod','linear'); 


addOptional(InPar,'Plot',true); 

parse(InPar,varargin{:});
InPar = InPar.Results;

if isempty(InPar.TargetList)
    % simulation mode
   
    InPar.DecRange = [-30 90]./RAD;
    
    [TileList,TileArea]=celestial.coo.tile_the_sky(56,42);
    Flag = TileList(:,2)>=InPar.DecRange(1) & TileList(:,2)<=InPar.DecRange(2);
    TileList = TileList(Flag,:);
    Ntarget  = size(TileList,1);
    %             RA, Dec, LastObs, MainCounter, NightCounter
    InPar.TargetList = [TileList(:,1:2), zeros(Ntarget,3)];
    
    InPar.JD = celestial.time.julday([1 3 2021 21 00 00]);
    
end


% get JD
if isempty(InPar.JD)
    % get JD of now
    JD = celestial.time.julday;
else
    JD = InPar.JD;
end

% TargetList
TargetList = InPar.TargetList;

% get visibility
if isempty(InPar.Visibility)
    [Res, LimitAM]     = celestial.scheduling.coo_visibility(JD,TargetList(:,1),TargetList(:,2),'AllNightVisibility',true,...
                                                            'Lon',InPar.Lon,'Lat',InPar.Lat,...
                                                            'TimeStep',InPar.TimeStep,...
                                                            'SunAltLimit',InPar.SunAltLimit,...
                                                            'MaxAM',InPar.MaxAM,...
                                                            'AzAltLimit',InPar.AzAltLimit,...
                                                            'MinMoonDistIllum',InPar.MinMoonDistIllum,...
                                                            'MinVisibilityTime',InPar.MinVisibilityTime,...
                                                            'FactorVisibilityAM',InPar.FactorVisibilityAM,...
                                                            'CleanTime',true,...
                                                            'InterpMethod',InPar.InterpMethod);
else
    Res        = InPar.Visibility;
    LimitAM    = InPar.MaxAM;
end


[W,ResW] = celestial.scheduling.weight_cadence(JD,TargetList(:,3),TargetList(:,5),'NightCadence',InPar.NightCadence,...
                                                                                  'MainCadence',InPar.MainCadence,...
                                                                                  'Nfast',InPar.Nfast,...
                                                                                  'MainWFun',InPar.MainWFun,...
                                                                                  'NightWFunPar',InPar.NightWFunPar,...
                                                                                  'NightWFun',InPar.NightWFun,...
                                                                                  'NightWFunPar',InPar.NightWFunPar);


W = W.';
%visibility window to end of night

% search for current time in Visibility matrix (Res)
[TimeDiff,Itime] = min(abs(Res.VecJD-JD));
if TimeDiff>InPar.TimeStep
    % Time is out of visibility matrix
    % error

    Status  = false;
    Problem = 'Time is out of visibility matrix';
else
    Status  = true;
    Problem = '';
end

% total time during the night to observe each target
TimeLeftForTarget = sum(Res.Flag(Itime:end,:),1).*InPar.TimeStep;  % [day]
TimeLeftForTargetLim = sum(Res.FlagLim(Itime:end,:),1).*InPar.TimeStep;  % [day]

% number of observations still needed for each target
NobsLeft          = InPar.Nfast - TargetList(:,5);

% Time needed to complete the target observations for the night
TimeNeedeForTarget = (NobsLeft-1).*InPar.NightCadence + (JD - TargetList(:,3));
% indicated targets with at least one obs during the night
FlagTargetStarted = TargetList(:,5)>0;

TimeNeedeForTarget(~FlagTargetStarted) = TimeLeftForTargetLim(~FlagTargetStarted);
TimeNeedeForTarget = TimeNeedeForTarget.';
FlagEnoughTime = TimeNeedeForTarget<TimeLeftForTarget & TimeLeftForTarget>0;



%WV = Res.Flag(Itime,:).*W.*FlagEnoughTime;
WV = Res.FlagLim(Itime,:).*W.*FlagEnoughTime;
Itarget = find(WV>0);

Nfound = numel(Itarget);
[~,SI]=sort(WV,'descend');
Ind = SI(1:Nfound);

% current visibility statistics of selected objects
CurVis = celestial.scheduling.coo_visibility(JD,TargetList(:,1),TargetList(:,2),'AllNightVisibility',false,...
                                                    'Lon',InPar.Lon,'Lat',InPar.Lat,...
                                                    'TimeStep',InPar.TimeStep,...
                                                    'SunAltLimit',InPar.SunAltLimit,...
                                                    'MaxAM',InPar.MaxAM,...
                                                    'AzAltLimit',InPar.AzAltLimit,...
                                                    'MinMoonDistIllum',InPar.MinMoonDistIllum,...
                                                    'MinVisibilityTime',InPar.MinVisibilityTime,...
                                                    'FactorVisibilityAM',InPar.FactorVisibilityAM,...
                                                    'CleanTime',true,...
                                                    'InterpMethod',InPar.InterpMethod);

% prepare output
ResS.Visibility = Res;
ResS.CurVis     = CurVis;
ResS.Nfound     = Nfound;
ResS.Ind        = Ind;



