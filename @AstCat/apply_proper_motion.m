function [OutCat]=apply_proper_motion(AstC,varargin)
%--------------------------------------------------------------------------
% apply_proper_motion function                               class/@AstCat
% Description: Apply proper motion parallax and radial velocity to
%              RA/Dec in an AstCat object.
% Input  : - An AstCat object that contains RA, Dec and optionally
%            PM_RA, PM_Dec, Plx, and RV columns.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'EpochInRA'  - Input epoch for RA. If empty then will assume
%                           the input is a SIM object and will attempt
%                           to use julday.m to get EpochIn.
%                           If this is a vector then each element will
%                           be applied to each element in the AstCat
%                           object.
%                           Alternatively this can be a string containing
%                           the column name of the EpochIn in RA.
%                           Default is empty.
%            'EpochInDec' - Input epoch for Dec. If empty then will use the
%                           values from EpochInRA.
%                           If this is a vector then each element will
%                           be applied to each element in the AstCat
%                           object.
%                           Alternatively this can be a string containing
%                           the column name of the EpochIn in Dec.
%                           Default is empty.
%            'EpochInUnits'-Units of EpochIn. Options are:
%                           'jd' - Julian days. Default.
%                           'jy' - Julian years.
%                           'by' - Bessilian years.
%            'EpochOut'   - Output epoch. Defualt is the cuurent time.
%            'EpochOutUnits'-Units of EpochOut. Options are:
%                           'jd' - Julian days. Default.
%                           'jy' - Julian years.
%                           'by' - Bessilian years.
%            'ColRA'      - Column name containing the RA information.
%                           If ColUnits field is empty then will assume
%                           RA is given in radians.
%                           Default is 'RA'.
%            'ColDec'     - Column name containing the Dec information.
%                           If ColUnits field is empty then will assume
%                           Dec is given in radians.
%                           Default is 'Dec'.
%            'ColPM_RA'   - Column name containing the PM_RA information.
%                           If ColUnits field is empty then will assume
%                           PM_RA is given in 'mas/yr'.
%                           Default is 'pmRA'.
%            'ColPM_Dec'  - Column name containing the PM_Dec information.
%                           If ColUnits field is empty then will assume
%                           PM_Dec is given in 'mas/yr'.
%                           Default is 'pmDec'.
%            'ColPlx'     - Column name containing the Plx information.
%                           If ColUnits field is empty then will assume
%                           Plx is given in 'mas'.
%                           Default is 'Plx'.
%            'ColRV'      - Column name containing the RV information.
%                           If ColUnits field is empty then will assume
%                           RV is given in 'km/s'.
%                           Default is 'RV'.
%            'NewColName' - Cell array of two new column names in which
%                           the RA/Dec will be stored. If empty, then
%                           will store RA/Dec in the old columns.
%                           Default is empty.
% Output : - An AstCat object that contains the RA/Dec at the requested
%            output epoch.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [OutCat]=apply_proper_motion(RefCat);
% Reliable: 2
%--------------------------------------------------------------------------


DefV.EpochInRA         = [];     % attempt to read from SIM header
DefV.EpochInDec        = [];     
DefV.EpochInUnits      = 'JD';   % JY/year/yr, BY, JD, 
DefV.EpochOut          = celestial.time.julday;
DefV.EpochOutUnits     = 'JD';   % JY/year/yr, BY, JD, 
DefV.ColRA             = 'RA';
DefV.ColDec            = 'Dec';
DefV.ColPM_RA          = 'pmRA';
DefV.ColPM_Dec         = 'pmDec';
DefV.ColPlx            = 'Plx';
DefV.ColRV             = 'RV';
%!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!
%logical parallax correction 
DefV.ApplyParallax     = false;
%!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!
DefV.NewColName        = [];   % if empty then insert instead of existing RA/Dec cols

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


NepochinRA  = numel(InPar.EpochInRA);
NepochinDec = numel(InPar.EpochInDec);
Nepochout   = numel(InPar.EpochOut);

Ncat   = numel(AstC);
OutCat = AstC;
for Icat=1:1:Ncat
    
    % EpochIn of RA
    if (isempty(InPar.EpochInRA))
        EpochInRA    = julday(AstC);
        EpochInUnits = 'JD';
    else
        [EpochInRA,EpochInUnits] = col_get(AstC(Icat),InPar.EpochInRA,InPar.EpochInUnits);
    end
    % convert to JD
    EpochInRA = celestial.time.year2jd(EpochInRA,InPar.EpochInUnits);
    
    % EpochIn of Dec
    if (isempty(InPar.EpochInDec))
        EpochInDec = EpochInRA;
    else
        [EpochInDec,~] = col_get(AstC(Icat),InPar.EpochInDec,InPar.EpochInUnits);
    end
    % convert to JD
    EpochInDec = celestial.time.year2jd(EpochInDec,InPar.EpochInUnits);
    
    
    % convert EpochOut from years to JD
    if (ischar(InPar.EpochOut))
        % EpochOut is a column in AstC
        ColEpochOut = colname2ind(AstC(Icat),InPar.EpochOut);
        EpochOut = AstC(Icat).Cat(:,ColEpochOut);
    else
        % EpochOut is a scalar or vector
        EpochOut = celestial.time.year2jd(InPar.EpochOut(Nepochout),InPar.EpochOutUnits);
    end
    
    % column indices [RA, Dec, PM_RA, PM_Dec, Plx, RV]
    ColInd = colname2ind(AstC(Icat),{InPar.ColRA, InPar.ColDec,...
                                     InPar.ColPM_RA, InPar.ColPM_Dec,...
                                     InPar.ColPlx, InPar.ColRV});
    
    if (any(isnan(ColInd)))
        % can not calculate PM
        % copy as is
        warning('Proper motion information is not in catalog');
    else
        
        Ncol = size(AstC(Icat).Cat,2);
        if (isempty(AstC(Icat).ColUnits))
            warning('Units are not provided in catalog - assume units are correct');
            % populate ColUnits with default units
            AstC(Icat).ColUnits = cell(1,Ncol);
            DefUnits = {'rad','rad','mas/yr','mas/yr','mas','km/s'};
            FlagNN = ~isnan(ColInd);
            AstC(Icat).ColUnits(ColInd(FlagNN)) = DefUnits(FlagNN); 
        end



        % RA/Dec
        if (isnan(ColInd(1)) || isnan(ColInd(2)))
            error('RA/Dec columns are not in catalog');
        end
        CI = 1;    % numeber of element in ColInd
        RA   = AstC(Icat).Cat(:,ColInd(CI));
        %RA  = RA.*convert_units(AstC(Icat).ColUnits{CI},'rad');
        RA  = RA.*convert.angular(AstC(Icat).ColUnits{CI},'rad');
        CI = 2;    % numeber of element in ColInd
        Dec  = AstC(Icat).Cat(:,ColInd(2));
        %Dec = Dec.*convert_units(AstC(Icat).ColUnits{CI},'rad');
        Dec = Dec.*convert.angular(AstC(Icat).ColUnits{CI},'rad');

        Nrow = numel(RA);

        % PM_RA
        CI = 3;   % numeber of element in ColInd
        if (isnan(ColInd(CI)))
            PM_RA = zeros(Nrow,1);
        else
            PM_RA = AstC(Icat).Cat(:,ColInd(CI));
            %PM_RA = PM_RA.*convert_units(AstC(Icat).ColUnits{ColInd(CI)},'mas/yr');
            PM_RA = PM_RA.*convert.units(AstC(Icat).ColUnits{ColInd(CI)},'mas/yr');
        end


        % PM_Dec
        CI = 4;   % numeber of element in ColInd
        if (isnan(ColInd(CI)))
            PM_Dec = zeros(Nrow,1);
        else
            PM_Dec = AstC(Icat).Cat(:,ColInd(CI));
            %PM_Dec = PM_Dec.*convert_units(AstC(Icat).ColUnits{ColInd(CI)},'mas/yr');
            PM_Dec = PM_Dec.*convert.units(AstC(Icat).ColUnits{ColInd(CI)},'mas/yr');
        end


        % Plx
        CI = 5;   % numeber of element in ColInd
        if (isnan(ColInd(CI)))
            Plx = zeros(Nrow,1);
        else
            Plx = AstC(Icat).Cat(:,ColInd(CI));
            %Plx = Plx.*convert_units(AstC(Icat).ColUnits{ColInd(CI)},'mas');
            Plx = Plx.*convert.angular(AstC(Icat).ColUnits{ColInd(CI)},'mas');
        end


        % RV
        CI = 6;   % numeber of element in ColInd
        if (isnan(ColInd(CI)))
            RV = zeros(Nrow,1);
        else
            RV = AstC(Icat).Cat(:,ColInd(CI));
            %RV = RV.*convert_units(AstC(Icat).ColUnits{ColInd(CI)},'km/s');
            RV = RV.*convert.units(AstC(Icat).ColUnits{ColInd(CI)},'km/s');
        end


        % calculate RA/Dec at EpochOut
        RV(isnan(RV)) = 0;
        % -------- !!!!!! ---------
        
        %Check for user flag for apply parallax barycentric
        if (InPar.ApplyParallax)
            [RA,Dec] = celestial.coo.proper_motion_parallax(EpochOut,EpochInRA,EpochInDec,RA,Dec,PM_RA,PM_Dec,Plx,RV);
        else
            [RA,Dec] = celestial.coo.proper_motion(EpochOut,EpochInRA,EpochInDec,RA,Dec,PM_RA,PM_Dec,Plx,RV);
        end
        % -------- !!!!!! ---------
        
        if (isempty(InPar.NewColName))
            % insert RA/Dec instead of existing columns
            OutCat(Icat) = col_replace(OutCat(Icat),[RA,Dec],ColInd(1:2));
        else
            % add new columns
            OutCat(Icat) = col_insert(OutCat(Icat),RA,Ncol+1);
            OutCat(Icat) = col_insert(OutCat(Icat),Dec,Ncol+2);
        end
    end
end