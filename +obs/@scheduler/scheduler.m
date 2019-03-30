%--------------------------------------------------------------------------
% AstCat class                                                       class
% Description: A class of structure array of catalogs (AstCat).
%              Note that this class is a subset of the SIM class.
% Input  : null
% Output : null
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef scheduler < handle
    properties (SetAccess = public)
        Targets        = AstCat(1,1); % AstCat of targets
        NextTarget                    % Index of next target
        FileName       = 'WFAST1_TargetList.mat';
        
    
    end
    
    properties (Constant = true)
        
    end
    
    properties (Hidden = true)
                     % Twilight -12, rise, set, Twilight -12
        SunMorningTwilight   = NaN;
        SunRiseSetLastUpdate = NaN;
        
        
    end
    
    properties (SetAccess = immutable)
        % property can be set only in the constructor.
        
        Lon       % observatory longitude [rad]
        Lat       % observatory latitute [rad]
        Height    % observatory height [m]
        Clock  
        
        SunPar      % Sun parameters for obs % MaxSunAltObs
                                             % MaxSunAltOpen 
        MoonPar     % Moon parameters for obs
        
        AltLimit    = 15;      % deg
        SunAltLimit = -11.5;   % deg
        
    end
    
    
    %-------------------
    %--- Constractor ---
    %-------------------
    methods
        
        function Sched=scheduler(varargin)
            % obs class constractor
            
            DefV.TargetsFile = [];
            DefV.Lon    = 0.60671659; % 34.7623;   
            DefV.Lat    = 0.53402537; % 30.5974; 
            DefV.Height = 875;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            Sched.Lon    = InPar.Lon;
            Sched.Lat    = InPar.Lat;
            Sched.Height = InPar.Height;
            
            % upload targets
            if (exist(InPar.TargetsFile,'file')>0)
                Sched.Targets = Util.IO.load2(InPar.TargetsFile);
            else
                % file doesn't exist
                AstC = obs.scheduler.create_targets_list;
                Sched.Targets = AstC;
            end
            
            
        end
        

    end
    
    methods (Static)
        function [Info,ColCell,ColUnits]=targets_colcell
            % Return default Targets AstCat ColCell (column names)
            Ind = 0;
            Ind = Ind + 1;
            Info(Ind).Col   = 'Name';
            Info(Ind).Units = '';
            Info(Ind).Def   = 'NoName';
            
            Ind = 0;
            Ind = Ind + 1;
            Info(Ind).Col   = 'Field';
            Info(Ind).Units = '';
            Info(Ind).Def   = NaN;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'RA';
            Info(Ind).Units = 'deg';
            Info(Ind).Def   = NaN;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'Dec';
            Info(Ind).Units = 'deg';
            Info(Ind).Def   = NaN;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'AddedTime';
            Info(Ind).Units = '';
            Info(Ind).Def   = NaN;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'IsToO';
            Info(Ind).Units = '';
            Info(Ind).Def   = false;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'Priority';
            Info(Ind).Units = '';
            Info(Ind).Def   = 1;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'Group';
            Info(Ind).Units = '';
            Info(Ind).Def   = 'WIS';
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'Project';
            Info(Ind).Units = '';
            Info(Ind).Def   = 'default';
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'ExpTime';
            Info(Ind).Units = 's';
            Info(Ind).Def   = 30;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'Nseq';  % number of continous exposures
            Info(Ind).Units = 's';
            Info(Ind).Def   = 1;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MaxNobs';
            Info(Ind).Units = '';
            Info(Ind).Def   = Inf;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'Nobs';
            Info(Ind).Units = '';
            Info(Ind).Def   = 0;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'LastTime';
            Info(Ind).Units = '';
            Info(Ind).Def   = NaN;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'StartTime';
            Info(Ind).Units = '';
            Info(Ind).Def   = 0;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'EndTime';
            Info(Ind).Units = '';
            Info(Ind).Def   = 0;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'NperNight';
            Info(Ind).Units = '';
            Info(Ind).Def   = 1;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'Cadence';
            Info(Ind).Units = '';
            Info(Ind).Def   = 10;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MaxPriority';
            Info(Ind).Units = '';
            Info(Ind).Def   = 0.9;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MinAM';
            Info(Ind).Units = '';
            Info(Ind).Def   = 1;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MaxAM';
            Info(Ind).Units = '';
            Info(Ind).Def   = 2.0;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MaxAbsHA';
            Info(Ind).Units = 'deg';
            Info(Ind).Def   = 90;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MinAlt';
            Info(Ind).Units = 'deg';
            Info(Ind).Def   = 15;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MaxAlt';
            Info(Ind).Units = 'deg';
            Info(Ind).Def   = 90;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MaxSunAlt';
            Info(Ind).Units = 'deg';
            Info(Ind).Def   = -11.5;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MaxMoonAlt';
            Info(Ind).Units = 'deg';
            Info(Ind).Def   = 90;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MinMoonDist';
            Info(Ind).Units = 'deg';
            Info(Ind).Def   = 20;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MinMoonIllum';
            Info(Ind).Units = '';
            Info(Ind).Def   = 0;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MaxMoonIllum';
            Info(Ind).Units = '';
            Info(Ind).Def   = 1;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MinMoonSkyExcess';
            Info(Ind).Units = 'deg';
            Info(Ind).Def   = 1;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MaxMoonSkyExcess';
            Info(Ind).Units = 'deg';
            Info(Ind).Def   = 1;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MinSeeing';
            Info(Ind).Units = 'arcsec';
            Info(Ind).Def   = 0;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MaxSeeing';
            Info(Ind).Units = 'arcsec';
            Info(Ind).Def   = 10;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MinTransp';
            Info(Ind).Units = '';
            Info(Ind).Def   = 0.1;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MaxTransp';
            Info(Ind).Units = '';
            Info(Ind).Def   = 1;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MinParAng';
            Info(Ind).Units = 'deg';
            Info(Ind).Def   = -Inf;
            
            Ind = Ind + 1;
            Info(Ind).Col   = 'MaxParAng';
            Info(Ind).Units = 'deg';
            Info(Ind).Def   = Inf;
            
            ColCell  = {Info.Col};
            ColUnits = {Info.Units};
          
        end
        
        function AstC=create_targets_list
            % Create an empty targets lits for scheduler object
            % Package: obs.scheduler
            % Description: Create an empty targets lits for scheduler object
            
            
            CatField      = AstCat.CatField;
            ColCellField  = AstCat.ColCellField;
            ColUnitsField = AstCat.ColUnitsField;
            
            AstC = AstCat(1,1);
            [~,ColCell,ColUnits] = obs.scheduler.targets_colcell;
            Ncol = numel(ColCell);
            
            AstC.(CatField) = array2table(ones(0,Ncol),'VariableNames',ColCell);
            AstC.(CatField).Properties.VariableUnits = ColUnits;
            AstC.(ColCellField)  = ColCell;
            AstC.(ColUnitsField) = ColUnits;
            AstC                 = colcell2col(AstC);
            
        end
        
        
        
    end
    
    methods
        function Sched=add_target(Sched,varargin)
            % add target to Scheduler object targets
            % Package: obs.scheduler
            % Description: add target to Scheduler object targets
            % Input  * Arbitrary number of pairs of ...,key,val,...
            %          arguments.
            %          Arguments names are listed in obs.scheduler.targets_colcell
            % Output : - Scheduler object with updated Target list
            
            CatField   = AstCat.CatField;
            
            [Info,ColCell] = obs.scheduler.targets_colcell;
            
            Ncol = numel(ColCell);
            for Icol=1:1:Ncol
                DefV.(Info(Icol).Col) = Info(Icol).Def;
            end
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);
            
            % add the values to the table
            Ntarget = size(Sched.Targets.(CatField),1);
            RowCell = cell(1,Ncol);
            for Icol=1:1:Ncol
                RowCell{Icol} = InPar.(Info(Icol).Col);
                %Sched.Targets.(CatField).(Info(Icol).Col)(Ntarget+1) = InPar.(Info(Icol).Col);
            end
            Sched.Targets.(CatField) = [Sched.Targets.(CatField); RowCell];
            % update targets file
            Targets = Sched.Targets;
            save(Sched.FileName,'Targets','-v7.3'); 
            
        end
        
    end
    
    methods
        function [FlagAll,Flag] = is_observable(Sched)
            % Select observable targets
            % Package: obs.scheduler
            % Description: Given a scheduler object containing a list of
            %              targets, test which target is observable based
            %              on various criteria as indicated by the
            %              scheduler.Targets list.
            % Input  : - Scheduler object
            % Output : - A vector of logical flags for each line in Targets
            %            list indicating of object is observable according
            %            to its observability criteria.
            %          - Structure of flags for individual observability
            %            criteria.
            
            
            RAD        = 180./pi;
            SEC_IN_DAY = 86400;
            CatField   = AstCat.CatField;
            
            Ntarget    = size(Sched.Targets.(CatField),1);
            
            [~,ColCell,ColUnits] = obs.scheduler.targets_colcell;
            Col = cell2struct(num2cell(1:1:numel(ColCell)),ColCell,2);
            
            % get current JD
            JD = celestial.time.julday;
            InfoSun  = celestial.SolarSys.get_sun(JD,[Sched.Lon, Sched.Lat]);
            InfoMoon = celestial.SolarSys.get_moon(JD,[Sched.Lon, Sched.Lat]);
            LST      = celestial.time.lst(JD,Sched.Lon); % frac of day
            
       
            if isnan(Sched.SunRiseSetLastUpdate)
                % update Sun rise/set time
                %RiseSet = celestial.SolarSys.sun_rise_set(JD,[Sched.Lon, Sched.Lat, Sched.Height],0);
                VecJD   = JD + (0:1:1440).'./1440;
                [SunRA,SunDec] = celestial.SolarSys.suncoo(VecJD,'a');
                SunAzAlt       = celestial.coo.horiz_coo([SunRA, SunDec], VecJD, [Sched.Lon, Sched.Lat],'h');
                SunCr          = Util.find.find_local_zeros(VecJD, SunAzAlt(:,2) - Sched.SunAltLimit./RAD);
                if (isempty(SunCr))
                    % no tweilight in 24 hr
                    Sched.SunMorningTwilight = NaN;
                else
                    SunCr = SunCr(SunCr(:,2)>0,:);
                    if (isempty(SunCr))
                        % no tweilight in 24 hr
                        Sched.SunMorningTwilight = NaN;
                    else
                        Sched.SunMorningTwilight = SunCr(1,1);
                    end
                end
                
            end
        
            if isnan(Sched.SunMorningTwilight)
                % no tweilight - assume >12hr
                TimeTillEndOfNight = 12.*3600;    % seconds
            else
                TimeTillEndOfNight = (JD - Sched.SunMorningTwilight).*SEC_IN_DAY;  % seconds
            end
            
            RA        = Sched.Targets.(CatField).RA  .* convert.angular(ColUnits{Col.RA},'rad');
            Dec       = Sched.Targets.(CatField).Dec .* convert.angular(ColUnits{Col.Dec},'rad');
            HA        = LST.*2.*pi - RA;   % [rad]
            ParAng    = celestial.coo.parallactic_angle([RA, Dec], LST, Sched.Lat);
            HorizCoo  = celestial.coo.horiz_coo([RA, Dec], JD, [Sched.Lon, Sched.Lat],'h');
            Az        = HorizCoo(:,1);
            Alt       = HorizCoo(:,2);
            AM        = celestial.coo.hardie(pi./2 - Alt);
            MoonDist  = celestial.coo.sphere_dist(RA,Dec,InfoMoon.RA,InfoMoon.Dec);
            MoonIllum = celestial.SolarSys.moon_illum(JD);
            
            Flag.AltLimit           = Alt > Sched.AltLimit;
            Flag.TimeTillEndOfNight = (Sched.Targets.(CatField).ExpTime .* Sched.Targets.(CatField).Nseq) < TimeTillEndOfNight;
            Flag.AM                 = AM >= Sched.Targets.(CatField).MinAM & ...
                                      AM <= Sched.Targets.(CatField).MaxAM;
            Flag.HA                 = abs(HA) < Sched.Targets.(CatField).MaxAbsHA;
            Flag.Nobs               = Sched.Targets.(CatField).Nobs <= Sched.Targets.(CatField).MaxNobs;
            Flag.Alt                = Alt > Sched.Targets.(CatField).MinAlt & ...
                                      Alt < Sched.Targets.(CatField).MaxAlt;
            Flag.SunAlt             = InfoSun.Alt < Sched.Targets.(CatField).MaxSunAlt;
            Flag.MoonAlt            = InfoMoon.Alt < Sched.Targets.(CatField).MaxMoonAlt;
            Flag.MoonDist           = MoonDist > Sched.Targets.(CatField).MinMoonDist;
            Flag.MoonIllum          = MoonIllum > Sched.Targets.(CatField).MinMoonIllum & ...
                                      MoonIllum < Sched.Targets.(CatField).MaxMoonIllum;
            Flag.ParAng             = ParAng > Sched.Targets.(CatField).MinParAng & ...
                                      ParAng < Sched.Targets.(CatField).MaxParAng;
     % NEED to update
            Flag.MoonSkyExcess      = true(Ntarget,1);
     % NEED to find a way to transmit seeing and transperancy
            %Flag.Seeing             
            %Flag.Transp
            
            % combine all observability flags
            FN      = fieldnames(Flag);
            Nfn     = numel(FN);
            FlagAll = true(Ntarget,1);
            for Ifn=1:1:Nfn
                FlagAll = FlagAll & Flag.(FN{Ifn});
            end
            
            
        end
        
        function [P] = calc_priority(Sched,varargin)
            % Calculate the cadence-based priority of all targets
            % Package: obs.scheduler
            % Description: Calculate the cadence-based priority of all
            %              targets.
            % Input  : - A Scheduler object with list of targets.
            % Output : - A vector of priority per target
            
            CatField = AstCat.CatField;
            
            DefV.CadenceFun           = 'bytime';
            DefV.AddPtoFirst          = 0.01; % for 'bytime' option - additional priority to add to first target in time window
            DefV.BasePriority         = 0.01;
            %DefV.TopPriority          = 1;
            DefV.ExtraPriority        = 0.1;
            DefV.RiseTime             = 0.1;  % In units of cadence
            DefV.TimeAtExtra          = 3;    % In units of cadence
            
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            t           = celestial.time.julday; 
            
            Cadence     = Sched.Targets.(CatField).Cadence;
            MaxPriority = Sched.Targets.(CatField).MaxPriority;
            switch lower(InPar.CadenceFun)
                case 'fermi'
                    %Cadence  = 
                    t        = celestial.time.julday; 
                    RiseTime = InPar.RiseTime.*Cadence;
                    P1 = (MaxPriority-InPar.ExtraPriority)./(1 + exp(-(t-Cadence)./RiseTime));
                    P2 = InPar.ExtraPriority.*exp(-(t-Cadence)./(Cadence.*InPar.TimeAtExtra));
                    P2(t<Cadence) = 0;
                    P  = P1 + P2;
                    
                case 'bytime'
                    % select targets by specific observing list
                    
                    Flag = t >= Sched.Targets.(CatField).StartTime & ...
                           t<Sched.Targets.(CatField).EndTime & ...
                           Sched.Targets.(CatField).LastTime < Sched.Targets.(CatField).StartTime & ...
                           Sched.Targets.(CatField).Nobs <= Sched.Targets.(CatField).MaxNobs;
                    P    = Sched.Targets.(CatField).MaxPriority .* Flag;
                    
                    % additional priority - prefer targets with earlier
                    % StartTime
                    tSt        = t - Sched.Targets.(CatField).StartTime;
                    tSt(tSt<0) = Inf;
                    [~,MinInd] = min(tSt);
                    
                    P(MinInd) = P(MinInd) + InPar.AddPtoFirst;
                    
                otherwise
                    error('Unknown CadenceFun option');
            end
            
            % check if target is observable
            P = P; %.*is_observable(Sched);
                    
            
            
        end
        
    end
    
    
    
end

            
