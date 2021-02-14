function [TimeHistory,History,ResS]=LAST_simulator(JD,varargin)
% Simulate LAST targets scheduling
% Package: +celestial.schedling
% Description:
% Input  : - JD of first night to simulate.
%          * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'DecRange' - Declination range of targets.
%                   Default is [-30 90].
%            'Ntel' - Number of telescopes to schedule. Default is 4.
%            'Nnight' - Number of nights to simulate. Default is 10.
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
%                   Default is true.
%            'MainCadence' - The main cadence of the survey [day]
%                       Default is 1.4.
%            'NightCadence' - The nightly cadence [day].
%                       Default is 40./1440.
%            'Nfast' - Number of epochs per night in the nightly cadence.
%                       Default is 2.
%            'MainWFun' - Weight function for main cadence.
%                       Default is @celestial.scheduling.fermiexp
%            'MainWFunPar' - Parameters for weight function for main
%                       cadence.
%                       Default is {1.4, 1, 0.03, 1, 0.5}.
%            'NightWFun' - Weight function for nightly cadence.
%                       Default is @celestial.scheduling.fermiexp
%            'NightWFunPar' - Parameters for weight function for nightly
%                       cadence.
%                       Default is {40./1440, 1, 0.003, 1.5, 0.5}.
%            'InterpMethod' - Tables interpolation. Default is 'linear'
%            'Plot' - Plot observed target as a function of time.
%                   Default is true.
% Output : - TimeHistory - A structure array with element per epoch,
%                   with the list of targets observed on each epoch.
%          - History - A structure arry with element per target.
%                   containing the target observing history.
%          - A structure with the followig fields:
%            .AllDiff - Vector of all time differences between targets,
%                   seperated with NaN between targets.
%            .AllAM - Vector of all observed AM.
%      By: Eran Ofek                      Oct 2020
% Example: [TimeHistory,History,ResS]=celestial.scheduling.LAST_simulator(celestial.time.julday([1 3 2021]),'Plot',false)
% histogram(ResS.AllAM,[1:0.05:2]')  % histogram of AM distribution
% plot number of vists per field
% axesm ('aitoff', 'Frame', 'on', 'Grid', 'on');
%scatterm(ResS.TargetList(:,2).*RAD,ResS.TargetList(:,1).*RAD,45,TL,'filled');
%H=colorbar;
%H.Label.String='Epochs/yr'; 
%H.Label.Interpreter='latex';
%G=celestial.coo.coco([(0:1:360)',zeros(361,1)],'g','j2000.0');
%hold on;
%plotm(G(:,2).*RAD,G(:,1).*RAD,'k.')
%axis off



RAD = 180./pi;

InPar = inputParser;
addOptional(InPar,'DecRange',[-30 90]./RAD);
%addOptional(InPar,'Ntel',8);
addOptional(InPar,'Ntel',4);

addOptional(InPar,'Nnight',10);
addOptional(InPar,'Lon',35./RAD);
addOptional(InPar,'Lat',30./RAD);

addOptional(InPar,'TimeStep',5./1440);  % day
addOptional(InPar,'SunAltLimit',-12./RAD);  
addOptional(InPar,'MaxAM',2);  
addOptional(InPar,'AzAltLimit',[0 15;90 15; 180 15; 270 15; 360 15]);  % [Az Alt] deg  
addOptional(InPar,'MinMoonDistIllum',[0 0; 0.1 1; 0.2 1; 0.3 1; 0.4 2; 0.5 3; 0.6 5;0.7 10;0.8 15; 0.9 30; 1.0 30]);  % [illum, MinDist]  

%addOptional(InPar,'MinVisibilityTime',5./24);  % [day] 
addOptional(InPar,'MinVisibilityTime',2./24);  % [day] 
addOptional(InPar,'FactorVisibilityAM',1.2);  % [day] 

%addOptional(InPar,'MainCadence',0.8);  % [day]
addOptional(InPar,'MainCadence',2.4);  % [day]
%addOptional(InPar,'NightCadence',30./1440); % [day]
addOptional(InPar,'NightCadence',40./1440); % [day]
%addOptional(InPar,'Nfast',2); % [day]
addOptional(InPar,'Nfast',8); % [day]

addOptional(InPar,'MainWFun',@celestial.scheduling.fermiexp); %@(t) 1.0+0.5.*exp(-t./1) ); % weight as a function of time since it is allowed to observe the target
%addOptional(InPar,'MainWFunPar', {0.8, 1, 0.03, 1, 0.5} );  %t0, DecayExp, SoftFermi, BaseW, ExtraW
addOptional(InPar,'MainWFunPar', {2.4, 1, 0.03, 1, 0.5} );  %t0, DecayExp, SoftFermi, BaseW, ExtraW
addOptional(InPar,'NightWFun',@celestial.scheduling.fermiexp); %@(t) 1.5+0.5.*exp(-t./1) ); % weight as a function of time since it is allowed to observe the target
%addOptional(InPar,'NightWFunPar',{30./1440, 1, 0.003, 1.5, 0.5});  %t0, DecayExp, SoftFermi, BaseW, ExtraW
addOptional(InPar,'NightWFunPar',{40./1440, 1, 0.003, 1.5, 0.5});  %t0, DecayExp, SoftFermi, BaseW, ExtraW

addOptional(InPar,'InterpMethod','linear'); 
addOptional(InPar,'Plot',true); 

parse(InPar,varargin{:});
InPar = InPar.Results;

[TileList,TileArea]=celestial.coo.tile_the_sky(56,42);
Flag = TileList(:,2)>=InPar.DecRange(1) & TileList(:,2)<=InPar.DecRange(2);
TileList = TileList(Flag,:);
Ntarget  = size(TileList,1);
%             RA, Dec, LastObs, MainCounter, NightCounter
TargetList = [TileList(:,1:2), zeros(Ntarget,3)];

TimeCounter = 0;
Hp = [];  % figure handle
NightJD = JD - 1;
for Inight=1:1:InPar.Nnight
    NightJD = NightJD + 1;
    [Res,LimitAM]     = celestial.scheduling.coo_visibility(NightJD,TileList(:,1),TileList(:,2),'AllNightVisibility',true,...
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

    Ntime   = numel(Res.VecJD);

    TargetList(:,5) = 0;
    for Itime=1:1:Ntime
        [Inight, Itime, Ntime]
        TimeCounter = TimeCounter + 1;
        
        JD = Res.VecJD(Itime);
        [W,ResW] = celestial.scheduling.weight_cadence(JD,TargetList(:,3),TargetList(:,5),'NightCadence',InPar.NightCadence,...
                                                                                          'MainCadence',InPar.MainCadence,...
                                                                                          'Nfast',InPar.Nfast,...
                                                                                          'MainWFun',InPar.MainWFun,...
                                                                                          'NightWFunPar',InPar.NightWFunPar,...
                                                                                          'NightWFun',InPar.NightWFun,...
                                                                                          'NightWFunPar',InPar.NightWFunPar);
                                                                                      
       

        W = W.';
        %visibility window to end of night
        
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
        if numel(Itarget)<InPar.Ntel
            %error('1')
            % Not enough targets - add targets
            numel(Itarget)
        end
        Nfound = min(numel(Itarget),InPar.Ntel);
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

        
        
        
        TargetList(Ind.',3) = JD;
        TargetList(Ind.',4:5) = TargetList(Ind.',4:5) + 1;
        
        TimeHistory(TimeCounter).JD = JD;
        TimeHistory(TimeCounter).RA = TargetList(Ind.',1).';
        TimeHistory(TimeCounter).Dec = TargetList(Ind.',2).';
        TimeHistory(TimeCounter).MainCounter = TargetList(Ind.',4).';
        TimeHistory(TimeCounter).NightCounter = TargetList(Ind.',5).';
        TimeHistory(TimeCounter).AM           = CurVis.AM(Ind);
        
        if TimeCounter>1
            CurRA  = TimeHistory(TimeCounter-1).RA;
            CurDec = TimeHistory(TimeCounter-1).Dec;
            NewRA  = TimeHistory(TimeCounter).RA;
            NewDec = TimeHistory(TimeCounter).Dec;
            
            %ResAs = celestial.scheduling.assign_targets2telescopes(NewRA,NewDec,CurRA,CurDec,'LST',Res.LST(Itime));
            
        end
        
        if InPar.Plot
            if ~isempty(Hp)
                Hp.Color = [0.5 0.5 0.5];
            else
                % first time
                axesm ('aitoff', 'Frame', 'on', 'Grid', 'on');
                hold on;
            end
            %Hp=plot(TimeHistory(TimeCounter).RA.*RAD,TimeHistory(TimeCounter).Dec.*RAD,'.','MarkerSize',14);
            %axis([0 360 -90 90]);
            if numel(Ind)>0
                Hp=plotm(TimeHistory(TimeCounter).Dec.*RAD,TimeHistory(TimeCounter).RA.*RAD,'.','MarkerSize',14);
            end
            
            drawnow;
            %M(TimeCounter) = getframe
            
        end
        
            
        
        for Ifound=1:1:Nfound
            History(Ind(Ifound)).JD(TargetList(Ind(Ifound),4)) = JD;
            History(Ind(Ifound)).AM(TargetList(Ind(Ifound),4)) = CurVis.AM(Ind(Ifound));
            History(Ind(Ifound)).MoonDist(TargetList(Ind(Ifound),4)) = CurVis.MoonDist(Ind(Ifound));
            
            
        end
        
        
        ResW;
    end
    Nodd(Inight) = sum(TargetList(:,5)./2~=floor(TargetList(:,5)./2));
    Neven(Inight) = sum(TargetList(:,5)./2==floor(TargetList(:,5)./2));
    
    for Itarget=1:1:Ntarget
        History(Itarget).NepochPerNight = TargetList(Itarget,5);
    end
    
    if InPar.Plot
        Hg = gca;
        LineObj = findobj(Hg.Children,'Type','Line');
        for Ip=1:1:numel(LineObj)
            LineObj(Ip).Color = [0.8 0.8 0.8]; 
            
        end
    end
    
end

AllDiff = [];
AllAM   = [];
AllMD   = [];
for I=1:1:numel(History)
    AllDiff = [AllDiff, NaN, diff(History(I).JD)];
    AllAM   = [AllAM, History(I).AM];
    AllMD   = [AllMD, History(I).MoonDist];
end

ResS.AllDiff = AllDiff;
ResS.AllAM   = AllAM;
ResS.TargetList = TargetList;




        
    