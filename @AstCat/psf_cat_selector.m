function [AstC,AllFlag]=psf_cat_selector(AstC,varargin)
%--------------------------------------------------------------------------
% psf_cat_selector function                                  class/@AstCat
% Description: Given a catalog of sources in image attempt to select
%              appropriate PSF stars based on some user defined criteria
%              including distance from image edge, peak flux, S/N, bit mask
%              values and more.
% Input  : - An AstCat object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'PosCol' - Source position columns in catalog.
%                       Default is {'XWIN_IMAGE','YWIN_IMAGE'}.
%            'BoundryDist' - Minimum allowed distance from image boundry
%                       (recomeded: use PSF half size).
%                       Default is 10 pix.
%            'SelectVal' - A three column cell array of
%                       {ColumnName, LowValue, HighValue}.
%                       For each column name (or arithmatic operation on
%                       column name) only entries between LowValue and
%                       HighValue will be selected.
%                       Default is {'PEAKF_VALTOT',-Inf,Inf;...
%                                   'SN',10,Inf}.
%            'SelectQuant' - A three column cell array of
%                       {ColumnName, LowQuantile, HighQuantile}.
%                       For each column name (or arithmatic operation on
%                       column name) only entries between the lower
%                       quantile and upper quantile of the column vales
%                       will be selected.
%                       Default is {'X2WIN_IMAGE',0.2,0.8;...
%                                   'Y2WIN_IMAGE',0.2,0.8;...
%                                   'XYWIN_IMAGE',0.2,0.8}.
%            'SatLevel' - Vector of saturation levels or an header keyword
%                       name or a cell array of header keyword names from
%                       which to obtain the saturation level.
%                       Default is 50000.
%            'ColMaxFlux' - Catalog column name containing the peak flux
%                       for the source. Default is 'FLUX_MAX'.
%            'ColSN'  - Catalog column name, or arithmatic operation over
%                       column names, that contain the source S/N.
%                       Default is 'FLUX_APER./FLUXERR_APER'.
%            'MinSN'  - Minimum S/N for a PSF star candidate.
%                       Default is 15.
%            'ColBitFlag' - Catalog column name containing the the bit mask
%                       flag for each source. Default is 'FLAGS'.
%            'MaskVal'- Bit mask value by which to reject bad sources.
%                       Default is 0 (i.e., no rejection).
%            'RejectWithNeigh' - Reject sources with neighbors.
%                       Default is true.
%            'FlagNeighPar' - Cell array of additional arguments to pass to
%                       AstCat/flag_neighbor. Default is {}.
% Output : - An AstCat object with the selected sources from the input
%            catalogs.
%          - A structure array of flags indicating which criteria where
%            satisfied by each source. The 'Flag' field is the combination of
%            all the criteria (i.e., the selection used to construct the
%            first output argument).
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [AstC,AllFlag]=psf_cat_selector(AstC);
% Reliable: 2
%--------------------------------------------------------------------------


DefV.PosCol             = {'XWIN_IMAGE','YWIN_IMAGE'};
DefV.BoundryDist        = 10;
DefV.SelectVal          = {'PEAKF_VALTOT',-Inf,Inf;...
                           'SN',10,Inf};
% DefV.SelectQuant        = {'X2WIN_IMAGE',0.2,0.8;...
%                            'Y2WIN_IMAGE',0.2,0.8;...
%                            'XYWIN_IMAGE',0.2,0.8};
DefV.SelectQuant        = {};

% DefV.SelectQuant        = {'X2WIN_IMAGE',0.2,1.8;...
%                            'Y2WIN_IMAGE',0.2,1.8;...
%                            'XYWIN_IMAGE',-0.2,0.2};                       
% DefV.SelectQuant        = {'X2WIN_IMAGE',0.0,0.3;...
%                            'Y2WIN_IMAGE',0.0,0.3;...
%                            'XYWIN_IMAGE',0.4,0.6};
%                        
DefV.SatLevel           = 60000;  % e- % or key name or cell of key names
DefV.ColMaxFlux         = 'PEAKF_VALTOT'; %'FLUX_MAX';
DefV.ColSN              = 'SN'; %FLUX_APER./FLUXERR_APER';   % 'SN'
DefV.MinSN              = 15;
DefV.ColBitFlag         = 'FLAGS';
DefV.MaskVal            = 0;    % reject if any of this bits are open - see bitmaks_find.m for options
DefV.RejectWithNeigh    = true;
DefV.FlagNeighPar       = {};

InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Nselval   = size(InPar.SelectVal,1);
Nselquant = size(InPar.SelectQuant,1);

Ncat    = numel(AstC);
Nrow    = sizecat(AstC);
AllFlag = Util.struct.struct_def({'Flag','Boundry','FluxMax','SN','Mask','Val','Quant','FlagNeigh'},Ncat,1);
for Icat=1:1:Ncat
    % for each catalog element
    
    
    %-----------------------------------------
    %--- remove stars near image boundries ---
    %-----------------------------------------
    if (~isempty(InPar.PosCol))
        CCDSEC  = ccdsec(AstC(Icat));
        CCDSECb = [CCDSEC(:,1) + InPar.BoundryDist + 1,...
                   CCDSEC(:,2) - InPar.BoundryDist - 1,...
                   CCDSEC(:,3) + InPar.BoundryDist + 1,...
                   CCDSEC(:,4) - InPar.BoundryDist - 1];
        [~,StFlag]=select_ccdsec(AstC(Icat),CCDSECb,InPar.PosCol);
        AllFlag(Icat).Boundry = StFlag.Flag;
    else
        % User did not request removal of stars near the boundry
        % select all stars
        AllFlag(Icat).Boundry = true(Nrow(Icat),1);
    end
    
    %-----------------------------------------------------
    %--- Remove stars using parameters value selection ---
    %-----------------------------------------------------
    AllFlag(Icat).Val = true(Nrow(Icat),1);
    if (~isempty(InPar.SelectVal))
        for Iselval=1:1:Nselval
            Val = col_arith(AstC(Icat),InPar.SelectVal{Iselval,1},'mat');
            AllFlag(Icat).Val = and(AllFlag(Icat).Val, Val>InPar.SelectVal{Iselval,2} & Val<InPar.SelectVal{Iselval,3});
        end
    end

    %--------------------------------------------------------
    %--- Remove stars using parameters quantile selection ---
    %--------------------------------------------------------
    AllFlag(Icat).Quant = true(Nrow(Icat),1);
    if (~isempty(InPar.SelectQuant))
        for Iselquant=1:1:Nselquant
            Val  = col_arith(AstC(Icat),InPar.SelectQuant{Iselquant,1},'mat');
            Low  = quantile(Val,InPar.SelectQuant{Iselquant,2});
            High = quantile(Val,InPar.SelectQuant{Iselquant,3});
            AllFlag(Icat).Quant = and(AllFlag(Icat).Quant, Val>Low & Val<High);
        end
    end


    %----------------------------------------------------
    %--- Remove staurated stars based on maximum flux ---
    %----------------------------------------------------
    if (~isempty(InPar.SatLevel))
        % Saturation level provided
        if any(strcmp(AstC(Icat).ColCell,InPar.ColMaxFlux))
            % MaxFlux column exist in catalog
            
            SatLevel = getkey_fromlist(AstC(Icat),InPar.SatLevel);
            SatLevel = SatLevel{1};
            if (~isnan(SatLevel))
                % Valid saturation level found
                
                FluxMax = col_get(AstC(Icat),InPar.ColMaxFlux);
                AllFlag(Icat).FluxMax = FluxMax<SatLevel;
            else
                warning('SatLevel not found - no MaxFlux selection');
                AllFlag(Icat).FluxMax = true(Nrow(Icat),1);
            end
        else
            warning('MaxFlux column was not found in catalog - no MaxFlux selection');
            AllFlag(Icat).FluxMax = true(Nrow(Icat),1);
        end
    else
        % User did not request removal of MaxFlux>SatLevel
        % select all stars
        AllFlag(Icat).FluxMax = true(Nrow(Icat),1);
    end
    
    %------------------------
    %--- Selection by S/N ---
    %------------------------
    if (~isempty(InPar.ColSN))
        % Select stars based on S/N
        SN = col_arith(AstC(Icat),InPar.ColSN,'mat',true);
        if (~all(isnan(SN)))
            AllFlag(Icat).SN   = SN>InPar.MinSN;
        else
            % all evaluated SN is NaN
            warning('All evaluated S/N are NaN - no S/N selection');
            AllFlag(Icat).SN = true(Nrow(Icat),1);
        end
    else
        % User did not request selection by S/N
        % select all stars
        AllFlag(Icat).SN = true(Nrow(Icat),1);
    end
    
    %-----------------------------
    %--- Selection by bit mask ---
    %-----------------------------
    if (InPar.MaskVal~=0)
        % select but bit mask
        if any(strcmp(AstC(Icat).ColCell,InPar.ColBitFlag))
            % ColBitFlag column exist in catalog
            AllFlag(Icat).Mask = bitand(col_get(AstC(Icat),InPar,ColBitFlag),InPar.MaskVal)==0;
        else
            warning('BitFlag column not in catalog - no bit mask selection');
            AllFlag(Icat).Mask = true(Nrow(Icat),1);
        end
        
    else
        % User did not request selection by bit mask
        % select all stars
        AllFlag(Icat).Mask = true(Nrow(Icat),1);
    end
    
    %-------------------------------
    %--- Selection by neighboors ---
    %-------------------------------
    % FFU
    if (InPar.RejectWithNeigh)
        AllFlag(Icat).FlagNeigh = ~flag_neighbor(AstC(Icat),InPar.FlagNeighPar{:});
    else
        AllFlag(Icat).FlagNeigh = true(Nrow(Icat),1);
    end
    
    %---------------------------------
    %--- Selection by point source ---
    %---------------------------------
    % FFU
    
    
    %---------------------------
    %--- Combining all flags ---
    %---------------------------
    AllFlag(Icat).Flag = AllFlag(Icat).Boundry & ...
                         AllFlag(Icat).FluxMax & ...
                         AllFlag(Icat).SN & ...
                         AllFlag(Icat).Mask & ...
                         AllFlag(Icat).Val & ...
                         AllFlag(Icat).Quant & ...
                         AllFlag(Icat).FlagNeigh;

    
    %--------------------------------------------------------------
    %--- Set the selected sources into the output AstCat object ---
    %--------------------------------------------------------------
    AstC(Icat) = col_select(AstC(Icat),Inf,AllFlag(Icat).Flag);
    
end

        
    
    
    
    
    
