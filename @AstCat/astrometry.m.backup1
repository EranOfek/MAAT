function [ResAst,OrigSim]=astrometry(Sim,varargin)
% Search astrometric solution to an image or catalog.
% Package: @AstCat
% Description: Search astrometric solution to an image or catalog.
%              Given a guess coordinates to an image or catalogs (of X/Y
%              positions), the image plate scale, and rough rotation,
%              the program attemt to identify the stars in the image
%              relative to stars in a reference catalog.
%              The matched stars are than used to solve for a
%              transformation between the X/Y coordinates and the RA/Dec.
%              The program is designed to be robust against failures. As
%              such the matching is doing in sub sections of the main image
%              so if some are failed there is still a chance to recover.
%              The fitting is done, by default, using orthogonal
%              polynomilas that makes the solution more stable.
% Input  : - An AstCat object or a SIM object.
%            An AstCat object that contains a catalog of sources in the
%            image, or a SIM object that contains an image and optionaly a
%            catalog of sources. If the catalog of sources in a SIM object
%            is not provided, then the sources will be extracted.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            '
% Output : - A structure array of astrometric solutions and quality
%            parameters for each image.
%          - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ResAst=astrometry(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;

CatField  = AstCat.CatField;
OrigColXY = {'dup_X','dup_Y'};



%--- Source extractor parameters ---
DefV.RePopCat           = false;
DefV.SrcExtractorProg   = @mextractor;
DefV.SrcExtractorPar    = {};
%--- Approximate guess position ---
DefV.Scale              = 1.01; % arcsec/pix
DefV.RA                 = [];   % rad,...
DefV.Dec                = [];   % rad,...
DefV.Equinox            = 2000.0;
DefV.KeyScale           = 'SCALE'; %
DefV.KeyRA              = {'RA','OBJRA','OBJRAD','CRVAL1'};
DefV.KeyDec             = {'DEC','OBJDEC','OBJDECD','CRVAL2'};
DefV.KeyEquinox         = {'EQUINOX'};
DefV.UnitsRA            = [];       % [], 'deg','r'
DefV.UnitsDec           = [];       % [], 'deg','r'
%--- Reference catalog ---
DefV.RefCat             = 'GAIADR2';  %@get_ucac4; %@wget_ucac4;   % string, function, struct
DefV.RCrad              = 0.8./RAD; %0.8/RAD;   % [radian]
DefV.RefCatMagRange     = [12 19.0]; %19.0];
DefV.UseMagRangeCat     = false;
DefV.MaxRefStars        = 15000;
DefV.Shape              = 'box';
DefV.RC_ColRA           = 'RA';
DefV.RC_ColDec          = 'Dec';
DefV.RC_ColMag          = 'Mag_G'; %'MagModel'; %'ModelMag';
DefV.RC_ColColor        = 'Mag_BP-Mag_RP';
DefV.ApplyPM            = true; %true; %true;
DefV.RC_ColPM_RA        = 'PMRA';
DefV.RC_ColPM_Dec       = 'PMDec';
DefV.RC_ColPlx          = 'Plx';
DefV.RC_ColRV           = 'RV';
DefV.RC_EpochInRA       = 'Epoch';
DefV.RC_EpochInDec      = 'Epoch';
DefV.RC_EpochInUnits    = 'yr';
DefV.CutRefCat          = true; %false;
%--- Catalog ---
DefV.ImSize             = [];  % [x,y]
DefV.ColXc              = {'XWIN_IMAGE','X','xpos'};
DefV.ColYc              = {'YWIN_IMAGE','Y','ypos'};
DefV.CatColMag          = 'MAG_PSF';  % Mag column name in SIM
%--- Cleaning ---
DefV.CleanLines         = true;
DefV.FlagLinesPar       = {};
DefV.CleanOverDense     = true;
DefV.FlagOverdensePar   = {};
%--- Matching ---
DefV.Flip               = [1 -1]; %1 1];
DefV.MatchRotMethod     = 'xcrot';     % 'scan' | 'xcrot'
DefV.HistDistEdges      = (12:3:300)';
DefV.MinRot             = -15;
DefV.MaxRot             = +15;
DefV.StepRot            = 0.5;
DefV.ReFindMatches      = true;
DefV.SearchRangeX       = [-1000 1000]; %[-1000 1000];
DefV.SearchRangeY       = [-1000 1000]; %[-1000 1000];
DefV.SearchRangeFactor  = 0.5;
DefV.SearchStepX        = 4; %3;
DefV.SearchStepY        = 4; %3;
DefV.SearchRad          = 4; %3;   % pix
DefV.MaxMethod          = 'adaptive'; %'sort';       % {'max1'|'sort'|'
DefV.MinNinPeak         = 5;            % min number of matches in peak - if <1 then fraction of ListCat
DefV.MaxPeaks           = 10;           % maximum number of peaks to return
DefV.MaxPeaksCheckI     = 5;     % 10      % out of the returned peaks - this is the maximum nuber of peaks to check
DefV.MaxPeaksCheckF     = 5;     % 10      % out of the returned peaks - this is the maximum nuber of peaks to check
DefV.SelectBest         = 'comb';    % {'N','std','meanerr','comb'}
DefV.FalseAlarm         = 1e-4;
DefV.RegionMaxConn      = 8;
DefV.AdjustStarDensity  = true;
DefV.MaxDensityAdjust   = 0.5; % 0.5;
DefV.NstarsRescaleBlock = 9000;     % if Nstars> then shrink blocksize by x2 and increase StepRot by x2
DefV.ReScaleFactor      = 2;
DefV.BlockSize          = [1024 1024]; %[512 512]; %[1024 1024]; %'full';
DefV.BufferSize         = 200;

%--- Fitting ---
%DefV.UseCase_TranC      = {'affine',             5}; %{'affine_tt_cheby2_4', 100; 'affine_tt_cheby2_3', 70; 'affine_tt',          20; 'affine',             5};
DefV.UseCase_TranC      = {'affine_tt_cheby2_4', 100; 'affine_tt_cheby2_3', 70; 'affine_tt',          10; 'affine',             5};
DefV.Niter              = 1;
DefV.SigClip            = 5;
DefV.MaxResid           = 1./3600;  % assuming coordinates in deg
%--- Analysis ---
DefV.AnalysisBlockSize  = [512 512];
%--- Verbose/plot ---
DefV.Verbose            = true;
DefV.Plot               = false;

%maximum Excess noise in the reference catalog (GAIA)
DefV.MaxExcessNoise     = 10;
%threshold for proper motion errors(GAIA)
DefV.MaxPMerr           = [];
%Logical for applying parallax correction (barycentric)
DefV.ApplyParallax      = false;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% number of images
Nsim = numel(Sim);

%------------------------
%--- source extractor ---
%------------------------
if (SIM.issim(Sim))
    IsCatPop = isfield_populated(Sim, CatField);
    if (~any(IsCatPop))
        Sim(~IsCatPop) = InPar.SrcExtractorProg(Sim(~IsCatPop),InPar.SrcExtractorPar{:});
    else
        if (InPar.RePopCat)
            Sim(~IsCatPop) = InPar.SrcExtractorProg(Sim(~IsCatPop),InPar.SrcExtractorPar{:});
        end
    end
else
    % assume SIM is an AstCat object with populated catalog
end

%-----------------------------------------------
%--- Convert user supplied coo user to J2000 ---
%-----------------------------------------------
if (~isempty(InPar.RA) && ~isempty(InPar.Dec))
    % User supplied RA/Dec
    RA  = celestial.coo.convertdms(InPar.RA,'gH','r');
    Dec = celestial.coo.convertdms(InPar.Dec,'gD','R');
    if (InPar.Equinox~=2000)
        Coo2000 = celestial.coo.coco([RA,Dec],sprintf('j%06.1f',InPar.Equinox),'j2000.0');  % Na'ama, 20180524
        RA  = Coo2000(:,1);
        Dec = Coo2000(:,2);
    end    
else
    % Read RA/Dec from header
    Out     = getcoo(Sim,'KeyRA',InPar.KeyRA,'KeyDec',InPar.KeyDec,'KeyEquinox',InPar.KeyEquinox,'OutUnits','rad');
    Coo2000 = celestial.coo.coco([[Out.RA].',[Out.Dec].'],sprintf('j%06.1f',Out(1).Equinox),'j2000.0');
    RA  = Coo2000(:,1);
    Dec = Coo2000(:,2);
end

%--- Read scale ---
if (isempty(InPar.Scale))
    Val = cell2mat(mgetkey(Sim,InPar.KeyScale));
    if any(isnan(Val))
        error('Some header scale keyword values are not available');
    end
    
    InPar.Scale = Val;
else
    InPar.Scale = InPar.Scale(:).*ones(Nsim,1);
end

% Debuging   
%Dec = Dec - 10./3600./RAD;
% [RA,Dec]=xy2coo(S,[1; 1024],[1; 2048]); [RA1,Dec1]=xy2coo(S1,[1; 1024],[1; 2048]); (RA-RA1).*RAD.*3600, (Dec-Dec1).*RAD.*360

%------------------
%--- Image size ---
%------------------
if (~isempty(InPar.ImSize))
    ImSize = InPar.ImSize;    % [X, Y]
else
    ImSize = imagesize(Sim);  % [X, Y]
end


OrigSim = Sim;

% astrometry for each image
FitRes     = Util.struct.struct_def({'Par','ParErr','Resid','ResidX','ResidY','ResidT',...
                         'Rrms','Chi2','Nobs','Npar','Ndof','rms','AstT','Flag','rmsAll','RrmsAll','Wrms','MagRef','Xref','Yref'},Nsim,1);
Log        = Util.struct.struct_def({'PatternMatch','FitMerge','FitFinal','Analysis'},Nsim,1);
ReturnFlag = false(Nsim,1);
for Isim=1:1:Nsim

    
    Scale = InPar.Scale(Isim);
    
    %!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!
    %add column with the order of the objects in the catalog 
    %Y location sorted in order find the objects that used in the fit
    
    
    %[~,ColXc]     = select_exist_colnames(Sim(Isim),InPar.ColXc(:));
    %[~,ColYc]     = select_exist_colnames(Sim(Isim),InPar.ColYc(:));
    
    %Sim(Isim).Cat(:,end+1)=(1:1:length(Sim(Isim).Cat(:,ColXc))).'; 

    Sim(Isim) = col_insert(Sim(Isim),...
                           (1:1:size(Sim(Isim).Cat,1)).',...
                           numel(Sim(Isim).ColCell)+1,...
                           'IndexSimYsorted');
    
    %Sim(Isim).Col.IndexSimYsorted=length(Sim(Isim).Cat(1,:)); 
    
    %Sim(Isim).ColCell{end+1}='IndexSimYsorted'; 
    
    
    %!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!
    %--------------------------
    %--- Reference Catalaog ---
    %--------------------------
    if (AstCat.isastcat(InPar.RefCat) || isstruct(InPar.RefCat))
        % RefCat was provided
        RefCat = InPar.RefCat;

        % convert to AstCat object
        RefCat = AstCat.struct2astcat(RefCat);
    else
        % External catalog was not provided
        % try to retrieve
        RefCat = VO.search.cat_cone(InPar.RefCat,RA(1),Dec(1),InPar.RCrad,'RadiusUnits','rad','OutType','astcat');
    end
    
    % what to do if RefCat is empty
    if (isempty(RefCat.Cat))
        % No reference catalog found - astrometric solution failed
        %ResAst(Isim) = [];
        warning('Reference catalog is empty - no solution found');
    else
        % RefCat is not empty
        
        if isa(InPar.RefCat,'function_handle')
            RefCat = InPar.RefCat;
        else
            % clean the GAIA catalog
            switch lower(InPar.RefCat)
                case 'gaiadr1'
                    % remove sources with excess noise >5 sigma and outside the
                    % mag range
                    InPar.RC_ColMag = 'MagG';   % override mag column
                    F = RefCat.(CatField)(:,8) < 80 & RefCat.(CatField)(:,5)> InPar.RefCatMagRange(1) & RefCat.(CatField)(:,5)< InPar.RefCatMagRange(2);
                    RefCat.(CatField) = RefCat.(CatField)(F,:);

                case 'gaiadr2'
                    MagG = col_get(RefCat,{InPar.RC_ColMag});
                    ExcessNoise = col_get(RefCat,{'ExcessNoise'});
                    F = ExcessNoise<InPar.MaxExcessNoise & MagG> InPar.RefCatMagRange(1) & MagG< InPar.RefCatMagRange(2);
                    RefCat.(CatField) = RefCat.(CatField)(F,:);


            end
        end
        
        %---------------------------
        %--- apply proper motion ---
        %---------------------------
        if (InPar.ApplyPM)
            % applay proper motion, RV and parallax to star positions
            EpochOut = julday(Sim);  % get JD of image - (image epoch)
            %%% ----- !!!!!!!!!!!  ----- !!!!!!!!!!!!
            
            %The field 'ApplyParallax' added to the apply proper motion
            %call
            RefCat = apply_proper_motion(RefCat,'EpochInRA',InPar.RC_EpochInRA,...
                                                'EpochInDec',InPar.RC_EpochInDec,...
                                                'EpochInUnits',InPar.RC_EpochInUnits,...
                                                'EpochOut',EpochOut,...
                                                'EpochOutUnits','JD',...
                                                'ColPM_RA',InPar.RC_ColPM_RA,...
                                                'ColPM_Dec',InPar.RC_ColPM_Dec,...
                                                'ColPlx',InPar.RC_ColPlx,...
                                                'ColRV',InPar.RC_ColRV, ...
                                                'ApplyParallax', InPar.ApplyParallax);
                                            
            %apply the limit on the PM error (if given by the user)
            if (~isempty(InPar.MaxPMerr))
                IndForPMerr= RefCat.Cat(:,9)<InPar.MaxPMerr;
                RefCat.Cat=RefCat.Cat(IndForPMerr,:);
            end
            %%% ----- !!!!!!!!!!!  ----- !!!!!!!!!!!!
        end

        % Generate a version of the reference catalog with only selected columns
        % RA, Dec, Mag, Color
        [RC_ColRA, RC_ColDec, RC_ColMag] = colname2ind(RefCat,{InPar.RC_ColRA, InPar.RC_ColDec, InPar.RC_ColMag});
        RC                = col_arith(RefCat,{InPar.RC_ColRA,InPar.RC_ColDec,...
                                                  InPar.RC_ColMag,InPar.RC_ColColor},...
                                      'astcat',true);
        RC.ColCell        = {InPar.RC_ColRA,InPar.RC_ColDec,'Mag','Color'};
        RC_ColMag         = 3;   % Mag in 3rd column
        RC                = colcell2col(RC);

        % applay MaxRefStars
        % If Ref catalogs have too many starts select only the brightest
        NrefStars = size(RC.(CatField),1);
        if (NrefStars>InPar.MaxRefStars)
            [~,IsM] = sort(RC.(CatField)(:,RC_ColMag));
            MaxMag  = RC.(CatField)(IsM(InPar.MaxRefStars),RC_ColMag);
            FlagMM  = RC.(CatField)(:,RC_ColMag)<MaxMag;
            RC.(CatField) = RC.(CatField)(FlagMM,:);
            
        end
            
        
        % calculate average density of sources in RefCat (RC):
        % this is used to adjust the number of stars in the catalog and reference.
        DensityRef = size(RC.Cat,1)./(pi.*(InPar.RCrad.*RAD).^2);   % stars/deg^2


        % project catalog coordinates from celestial sphere to plane
        % X,Y in units of pixels (using first guess pixel scale).
        % Note that the X=0,Y=0 position is the position defined
        % by the guess RA/Dec (RA, Dec parameters).
        % This are the "projection-plan coordinates"
        % relative to the native long and lat of the fiducial point
        [X,Y]=projection(RC,'tan',[RC_ColRA RC_ColDec],[RAD.*3600./Scale RA Dec],'rad');
        % add the projected X/Y as the 1st and 2nd columns in the reference catalog
        ColRef   = {'X','Y'};
        RC = col_insert(RC,X,1,'X');
        RC = col_insert(RC,Y,2,'Y');
        RC_ColMag         = 5;  % Mag is now in the 5th column!
        % sort the reference catalog by the Y position
        RC = sortrows(RC,'Y');
        % plot(RC.Cat(:,1),RC.Cat(:,2),'.')

        %------------------
        %--- Prep image ---
        %------------------
        CatColMagInd = colname2ind(Sim(Isim),InPar.CatColMag);
        
        % mag range
        if (InPar.UseMagRangeCat)
            FlagM = Sim(Isim).(CatField)(:,CatColMagInd)>InPar.RefCatMagRange(1) & ...
                    Sim(Isim).(CatField)(:,CatColMagInd)<InPar.RefCatMagRange(2);
            Sim(Isim).(CatField) = Sim(Isim).(CatField)(FlagM,:);
            
        end
        
        % clean sources out of image boundries
        [~,ColXc]     = select_exist_colnames(Sim(Isim),InPar.ColXc(:));
        [~,ColYc]     = select_exist_colnames(Sim(Isim),InPar.ColYc(:));
        
        FlagIn = Sim(Isim).(CatField)(:,ColXc)>0 & ...
                 Sim(Isim).(CatField)(:,ColXc)<ImSize(1) & ...
                 Sim(Isim).(CatField)(:,ColYc)>0 & ...
                 Sim(Isim).(CatField)(:,ColYc)<ImSize(2);
        Sim(Isim).(CatField) = Sim(Isim).(CatField)(FlagIn,:);
             
        % clean lines/rows
        % search for excess of sources on the same line/row and remove
        if (InPar.CleanLines)
            FlagCat = flag_bleeding_src(Sim(Isim),InPar.FlagLinesPar{:});
            Sim(Isim).(CatField) = Sim(Isim).(CatField)(~FlagCat.Flag,:);
        end

        % clean overdensity regions
        % search for excess f sources in some regions and remove
        if (InPar.CleanOverDense)
            FlagCat = flag_overdense_src(Sim,InPar.FlagOverdensePar{:});
            Sim(Isim).(CatField) = Sim(Isim).(CatField)(~FlagCat.Flag,:);
        end
        
        % adjust number of stars in Sim such that it will be roughly equal to
        % the star density in the reference image
        
        
        % Note that density adjustment can be done only if stars magnitudes
        % are provided
        % Check if magnitude is ok
        if (~isempty(CatColMagInd))
            if (numel(unique(Sim(Isim).(CatField)(:,CatColMagInd)))>1)
                ValidMag = true;
            else
                ValidMag = false;
            end
        else
            ValidMag = false;
        end
        if (~ValidMag && InPar.Verbose)
            warning('Invalid magnitudes - avoid using');
        end

        if (InPar.AdjustStarDensity && ValidMag)
            DensityCat = size(Sim(Isim).Cat,1)./prod(ImSize(Isim,:).*Scale./3600);
            
            DensityRatio = DensityRef./DensityCat;
            
            if (DensityRatio<1)
                
                % adjust the number of sources in the Sim image
                Fraction = max(DensityRatio, InPar.MaxDensityAdjust);
                Quant    = quantile(Sim(Isim).Cat(:,CatColMagInd),Fraction);
                Flag     = Sim(Isim).Cat(:,CatColMagInd)<Quant;
                Sim(Isim).Cat = Sim(Isim).Cat(Flag,:);
                
                
            else
                DensityRatio = 1./DensityRatio;
               
                % adjust the number of sources in the Ref 
                Fraction = max(DensityRatio, InPar.MaxDensityAdjust);
                Quant    = quantile(RC.Cat(:,RC_ColMag),Fraction);
                Flag     = RC.Cat(:,RC_ColMag)<Quant;
                RC.Cat   = RC.Cat(Flag,:);
            end
                
        end

        
        
        
        % define VecRot
        % define BlockSize and step size based on number of stars
        Nsrc   = size(Sim(Isim).(CatField),1);
        if (Nsrc>InPar.NstarsRescaleBlock)
            if (InPar.Verbose)
                fprintf('ReScale BlockSize and StepRot by %f\n',InPar.ReScaleFactor);
            end
            BlockSize = InPar.BlockSize./InPar.ReScaleFactor;
            VecRot = (InPar.MinRot:InPar.ReScaleFactor.*InPar.StepRot:InPar.MaxRot).';
        else
            BlockSize = InPar.BlockSize;
            VecRot = (InPar.MinRot:InPar.StepRot:InPar.MaxRot).';
        end

        % get X and Y column indices in SIM
        %[ColXc,ColYc] = colname2ind(Sim(Isim),{InPar.ColXc,InPar.ColYc});   
        %[~,ColXc]     = select_exist_colnames(Sim(Isim),InPar.ColXc(:));
        %[~,ColYc]     = select_exist_colnames(Sim(Isim),InPar.ColYc(:));
        

        % Construct a table with [X,Y] columns only

        SimCat = AstCat.sim2astcat(Sim(Isim));

        % Duplicate X/Y columns in SimCat
        % Original X/Y are named: 'dup_X' and 'dup_Y'
        % transformed X/Y are '*WIN_IMAGE'...
        SimCat     = col_duplicate(SimCat,[ColXc,ColYc],OrigColXY);

        % select sources in sub image 
        % If the solution is done in sub image then this generate
        % a catalog for each sub image (block).
        
        SubCat     = subcat_regional(SimCat,ImSize(Isim,:),BlockSize,InPar.BufferSize,[ColXc,ColYc]);
        Nsub       = numel(SubCat);
      
        ResBest    = Util.struct.struct_def({'MaxHistMatch','MaxHistShiftX','MaxHistShiftY',...
                            'ShiftX','ShiftY','Tran','Nmatch','IndRef','IndCat',...
                            'MatchedCat','MatchedRef','MatchedResid','StdResid',...
                            'Std','MeanErr','BestRot','BestFlip'},Nsub,1);
                        
       ShiftRes = nan(Nsub,2);
       ImCenter = ImSize(Isim,:).*0.5;
       for Isub=1:1:Nsub
           %Isub
            % for each sub region
        
            % match ref catalog with image catalog
            % find also the correct image flip and rotation.
            % The origin of the reference catalog (RC) is [RA, Dec].
            % The origin of the image is its geometric center.
            % Note that the output MatchedCat
            % contains the columns: {'XWIN_IMAGE'  'YWIN_IMAGE'  'dup_X'  'dup_Y'}
            
            
            % SubCat contains sources in a sub image region
            % next line transform its coordinates from image corner to
            % image center (i.e., similar to that of the reference catalog)
            SubCat(Isub).(CatField)(:,[ColXc, ColYc]) = SubCat(Isub).(CatField)(:,[ColXc, ColYc]) - ImCenter;  % +[50 50];
            
            switch lower(InPar.MatchRotMethod)
                case 'xcrot'
                    % cross match the rotation using histogram of
                    % distances/angles
                    
                    %Add a condition that the SubCat isn't empty, In the
                    %case of empty SubCat - skip
                    if any(any(SubCat(Isub).(CatField)))
                        ResRot = ImUtil.pattern.match_pattern_rot(SubCat(Isub).(CatField),RC.(CatField),...
                                                              'CatColX',ColXc,...
                                                              'CatColY',ColYc,...
                                                              'HistDistEdges',InPar.HistDistEdges,...
                                                              'CutRefCat',InPar.CutRefCat,...
                                                              'SearchRangeX',InPar.SearchRangeX,...
                                                              'SearchRangeY',InPar.SearchRangeY);
                    else
                        continue;
                    end
                    % go over all rotational possibilities
                    Nrot = size(ResRot.MatchedRot,1);
                    K = 0;
                    VecRotSel = [];
                    for Irot=1:1:Nrot
                        if (any(all(ResRot.MatchedRot(Irot,3:4)==InPar.Flip,2)))
                            K = K + 1;
                            VecRotSel(K) = ResRot.MatchedRot(Irot,1);
                        end
                    end
                    
                    if (isempty(VecRotSel))
                        % No rotation candidate solution found
                        % go back to scanning
                        VecRotSel = VecRot;
                    else
                    
                        Fvrs = VecRotSel<InPar.MaxRot | (VecRotSel>InPar.MinRot | VecRotSel>(360+InPar.MinRot));
                        VecRotSel = VecRotSel(Fvrs);
                        %ImUtil.pattern.match_pattern_rot finds minus the rotation
                        VecRotSel = -VecRotSel;
                    end
                case 'scan'
                    % match rotation by scanning all possible rotations
                    VecRotSel = VecRot;
                otherwise
                    error('Unknown MatchRotMethod option');
            end
           
            [Res,IndBest,H] = ImUtil.pattern.match_pattern_shift_rot(SubCat(Isub).(CatField),RC.(CatField),...
                                                        VecRotSel,...
                                                        'ColXc',ColXc,...
                                                        'ColYc',ColYc,...
                                                        'Flip',InPar.Flip,...
                                                        'CutRefCat',InPar.CutRefCat,...
                                                        'SearchRangeX',InPar.SearchRangeX,...
                                                        'SearchRangeY',InPar.SearchRangeY,...
                                                        'SearchRangeFactor',InPar.SearchRangeFactor,...
                                                        'SearchStepX',InPar.SearchStepX,...
                                                        'SearchStepY',InPar.SearchStepY,...
                                                        'Radius',InPar.SearchRad);
                                                    
                                             
            if (~isempty(Res) && ~isempty(IndBest))      
                ShiftRes(Isub,:) = [Res(IndBest).ShiftX, Res(IndBest).ShiftY];
                
                ResSub(Isub) = Res(IndBest);   
            end
            
       end
       
       % DEBUG
       %for Isub=1:1:8; plot(SubCat(Isub).Cat(:,2),SubCat(Isub).Cat(:,3),'.'); hold on; end
       %for Isub=1:1:8; plot(ResSub(Isub).MatchedCat(:,1),ResSub(Isub).MatchedCat(:,2),'.'); hold on; end
       
       MedianShift = nanmedian(ShiftRes,1);
       FlagShift   = abs(ShiftRes - MedianShift)< (max(InPar.SearchStepX,InPar.SearchStepY).*5);
       %FlagShift   = abs(ShiftRes - MedianShift)< (max(InPar.SearchStepX,InPar.SearchStepY).*50); % for LFC
       
       SubGood     = and(FlagShift(:,1),FlagShift(:,2));
       MatchedCat = [];
       MatchedRef = [];
       
       for Isub=1:1:Nsub
            % join the matched catalogs from the sub images
            if (SubGood(Isub))
                if (isempty(MatchedCat))
                    MatchedCat = ResSub(Isub).MatchedCat;
                    MatchedRef = ResSub(Isub).MatchedRef;
                else
                    MatchedCat = [MatchedCat; ResSub(Isub).MatchedCat];
                    MatchedRef = [MatchedRef; ResSub(Isub).MatchedRef];
                end
            end
       end
       
       
       if (isempty(MatchedRef) || isempty(MatchedCat))
           ResAst(Isim).ShiftRes = ShiftRes;
           ResAst(Isim).SubGood = any(SubGood);
       else
       
           Nmatch = size(MatchedRef,1);
           
           % select transformation based on number of sources
           Iuse = find([InPar.UseCase_TranC{:,2}]<Nmatch,1);
           if (isempty(Iuse))
               error('Number of stars (Nmatch=%d) is too low for solution',Nmatch);
           end
           TranC = InPar.UseCase_TranC{Iuse,1};

           % A factor for normalizing the X/Y coordinates to unity
           NormXY = max(ImSize(Isim,:));
           % Fit transformation
%            ResAst(Isim) = ImUtil.pattern.fit_transform(MatchedRef,MatchedCat,TranC,'ImSize',ImSize(Isim,:),...
%                                                                           'BlockSize',InPar.AnalysisBlockSize,...
%                                                                           'NormXY',NormXY,...
%                                                                           'PolyMagDeg',3,...
%                                                                           'StepMag',0.1,...
%                                                                           'Niter',InPar.Niter,...
%                                                                           'SigClip',InPar.SigClip,...
%                                                                           'MaxResid',InPar.MaxResid,...
%                                                                           'Plot',InPar.Plot);
                                                                      
         
           % Special treatment for PV distortions
           % convert coordinates from pixels to deg
           CD = [1 0; 0 1].*Scale./3600;
           %!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!
           %clear many appearence of the same object.
           [indexes,ia,ic]=unique(MatchedCat(:,Sim(Isim).Col.IndexSimYsorted));
           MatchedRef=MatchedRef(ia,:);
           MatchedCat=MatchedCat(ia,:);
           %!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!

           MatchedRefCD        = MatchedRef;
           MatchedRefCD(:,1:2) = [CD*MatchedRefCD(:,1:2)']';
           MatchedCatCD        = MatchedCat;
           MatchedCatCD(:,[ColXc, ColYc]) = [CD*MatchedCatCD(:,[ColXc, ColYc])']';
           
           ResAst(Isim) = ImUtil.pattern.fit_transform(MatchedRefCD,MatchedCatCD,TranC,'ImSize',ImSize(Isim,:),...
                                                                          'BlockSize',InPar.AnalysisBlockSize,...
                                                                          'PixScale',InPar.Scale,...
                                                                          'CooUnits','deg',...
                                                                          'NormXY',1,...
                                                                          'ColCatX',ColXc,...
                                                                          'ColcatY',ColYc,...
                                                                          'PolyMagDeg',3,...
                                                                          'StepMag',0.1,...
                                                                          'Niter',InPar.Niter,...
                                                                          'SigClip',InPar.SigClip,...
                                                                          'MaxResid',InPar.MaxResid,...
                                                                          'Plot',InPar.Plot);
           
           
           ResAst(Isim).ShiftRes   = ShiftRes;
           ResAst(Isim).SubGood    = any(SubGood);
           %!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!
           %add the catalog data of the used objects with the cols data
           TempAstCat=AstCat;
           TempAstCat.Cat = MatchedCat(ResAst(Isim).FlagMag,:);
           TempAstCat.Col=SimCat.Col; 
           TempAstCat.ColCell=SimCat.ColCell; 
           ResAst(Isim).AstCat= TempAstCat;
           %indexes vector of the used objects in the original catalog
           ResAst(Isim).IndexInSim1=unique(MatchedCat(ResAst(Isim).FlagMag,ResAst(Isim).AstCat.Col.IndexSimYsorted));
           ResAst(Isim).IndexInSimN= ResAst(Isim).IndexInSim1(ResAst.FlagG);
           ResAst(Isim).FlagMag=[];
           %!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!
           
           %---------------------------------------------------------
           %--- Convert the transformation to WCS header keywords ---
           %---------------------------------------------------------
           % ImCenter - Coordinate zero in the image X/Y coordinate system
           % NormXY   - Coordinate normalization
           % RA,Dec   - assumed center [rad]
           
           %W = ClassWCS.tranclass2wcs(ResAst(Isim).TranC,'CooCenter',[RA,Dec], 'ImCenter',ImCenter, 'NormXY',NormXY, 'Scale',Scale);
           if (nargout>1)
               W = ClassWCS.tranclass2wcs_tpv(ResAst(Isim).TranC,'CooCenter',[RA,Dec], 'ImCenter',ImCenter, 'NormXY',NormXY, 'Scale',Scale,'CD',CD);
               OrigSim(Isim) = wcs2head(W,OrigSim(Isim));
             
             %add WCS field
               ResAst(Isim).WCS=OrigSim(Isim).WCS;
               
               % Vancky
               %add or update WCS field in OriginSim
               % it seems WCS in SIM is inherited from superclass WorldCooSys
               % thus xy2coo for SIM call function in WorldCooSys, need to
               % fix? now we have to W = ClassWCS.populate(OrigSim); and call
               % xy2coo(W,[X,Y]);
               ResAst(Isim).WCS  = W;
               OrigSim(Isim).WCS = W;
             
           end
           
           
     
           %Res(Isim).plot_resmag = @(Res) semilogy(Res(
       end
        
    end
end
