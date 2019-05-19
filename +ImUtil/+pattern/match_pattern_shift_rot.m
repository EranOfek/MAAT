function [Res,IndBest,H2]=match_pattern_shift_rot(Cat,Ref,Rotation,varargin)
% Match X/Y coordinates in two shifted and rotated lists.
% Description: ImUtil.pattern
% Input  : - Catalog list, containing at least two columns of [X,Y].
%            Should be sorted by Y.
%            If this is an AstCat object then AstCat/pattern_match_shift
%            will be called.
%          - Reference list, containing at least two columns of [X,Y].
%            Sorted by Y.
%          - Vector of rotations to check [deg].
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ColXc' - Column name or column index of the X coordinate in
%                      the catalog. Default is 1.
%            'ColYc' - Column name or column index of the Y coordinate in
%                      the catalog. Default is 2.
%            'ColXr' - Column name or column index of the X coordinate in
%                      the reference. Default is 1.
%            'ColYr' - Column name or column index of the Y coordinate in
%                      the reference. Default is 2.
%            'FunCatSelect' - Function handle for selecting good sources in
%                      input Cat. Flag=fun(Cat,additional_arg), where Flag
%                      is a logical vector indicating good sources.
%                      Default is [].
%            'FunCatSelectPar' - Additional parameters to pass to the
%                      FunCatSelect function. Default is {}.
%            'FunRefSelect' - Function handle for selecting good sources in
%                      input Ref. Flag=fun(Ref,additional_arg), where Flag
%                      is a logical vector indicating good sources.
%                      Default is [].
%            'FunRefSelectPar' - Additional parameters to pass to the
%                      FunRefSelect function. Default is {}.
%            'SearchRangeX' - X shift search range (in coordinate units).
%                      If empty then set to half the max X coordinates in
%                      the input Cat.
%                      Default is [].       
%            'SearchRangeY' - Y shift search range (in coordinate units).
%                      If empty then set to half the max X coordinates in
%                      the input Cat.
%                      Default is [].
%            'SearchStepX' - X shift search step size (in coordinate units).
%                      Default is 2.
%            'SearchStepY' - Y shift search step size (in coordinate units).
%                      Default is 2.
%            'MaxMethod' - Method by which to use the best match
%                      solution/s. Options are:
%                      'sort' - sort by number of matches and return the
%                               best solutions. Default.
%                      'max1' - Return the soltion with the largest
%                               number of matches.
%            'MinNinPeak' - Minimum number of histogram matching for a
%                      valid solution. Default is 5.
%                      If value is <1, then this is the fraction of matches
%                      out of the total number of sources in the catalog.
%            'MaxPeaks' - Maximum number of histogram peaks (matches) to
%                      return. Default is 10.
%            'MaxPeaksCheck' - Out of the returned peaks, this is the
%                      maximum nuber of peaks to check. Default is 10.
%                      Checking includes trying to use the transformation
%                      to match all the sources.
%                      If 0 then 'SelectBest' will use the MaxHistMatch
%                      value.
%            'CooType' - Coordinate type {'plane'|'sphere'}.
%                      Default is 'plane'.
%            'SearchRad' - Search radius for final matching. Default is 2.
%            'SearchMethod' - Catalg matching search method.
%                      See search_cat.m for options. Default is 'binms'.
%            'SelectBest' - How to select the best solution. Options are:
%                      'meanerr' - Choose the solution with the minimal
%                                  error on the mean.
%                      'N'       - Choose the solution with the largest
%                                  number of matched sources.
%                      'std'     - Choose the solution with the minmum
%                                  std.
%                      'comb'    - Use both meanerr and N. Default.
%            'FalseAlarm' - In the adaptive peak detection option,
%                      this is the false alarm probability that will be
%                      used in order to calculate the detection threshold.
%                      Default is 1e-4.
%            'RegionMaxConn' - connectivity for imregionalmax. Default is 8.
% Output : - Structure array of all the possible shift solutions.
%            The following fields are available:
%          - Index of the best shift solution in the structure of all
%            possible solutions.
%          - The 2D histogram of all possible combination of X and Y
%            distances.         
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res,IndBest,H2]=ImUtil.pattern.match_pattern_shift_rot
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;

if (nargin==0)
    % simulation mode
    Rotation = (-1:0.1:1)';
    Nstar = 1000;
    Ref = rand(Nstar,2).*2048;
    Ref = sortrows(Ref,2);
    Noverlap = 50;
    Cat = [Ref(1:Noverlap,1), Ref(1:Noverlap,2)];
    Cat = [Cat; rand(Nstar-Noverlap,2).*2048];
    Cat(:,1) = Cat(:,1) + 20 + randn(Nstar,1).*0.3;
    Cat(:,2) = Cat(:,2) + 30 + randn(Nstar,1).*0.3;
    Cat      = sortrows(Cat,2);
end



DefV.Radius           = 2;
DefV.Flip             = [1 1];
DefV.ColXc            = 1;            % or string
DefV.ColYc            = 2;
DefV.ColXr            = 1;
DefV.ColYr            = 2;
DefV.FunCatSelect     = [];
DefV.FunCatSelectPar  = {};
DefV.FunRefSelect     = [];
DefV.FunRefSelectPar  = {};
DefV.CutRefCat        = true;
DefV.SearchRangeX     = [-1000 1000];
DefV.SearchRangeY     = [-1000 1000];
DefV.SearchRangeFactor= 0.5;
DefV.SearchStepX      = 3;
DefV.SearchStepY      = 3;
DefV.MaxMethod        = 'adaptive'; %'sort';       % {'max1'|'sort'|'
DefV.FalseAlarm       = 1e-2;
DefV.MinNinPeak       = 5;            % min number of matches in peak - if <1 then fraction of ListCat
DefV.RetSolFracMax    = 0.5;          % return only solutions that their HistNmatch is larger than this fraction*max(HistNmatch)
DefV.MaxPeaks         = 10;            % maximum number of peaks to return in individual rotations
DefV.MaxPeaksCheck    = 10;            % out of the returned peaks - this is the maximum nuber of peaks to check
%DefV.SearchRad        = 2;
DefV.SortResBy        = 'N';
DefV.SelectBest       = 'comb';    % {'N','std','meanerr','comb'}
%DefV.RegionMaxConn = 8;

if (~isempty(varargin))
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);
else
    InPar = DefV;
end


% set SearchRangeX/Y
if (isempty(InPar.SearchRangeX))
    SizeX = max(Cat(:,InPar.ColXc)).*InPar.SearchRangeFactor;
    InPar.SearchRangeX = [-SizeX +SizeX];
end
if (isempty(InPar.SearchRangeY))
    SizeY = max(Cat(:,InPar.ColYc)).*InPar.SearchRangeFactor;
    InPar.SearchRangeY = [-SizeY +SizeY];
end

% If SearchRange is a scalar then assume the range is symmetric
if (numel(InPar.SearchRangeX)==1)
    InPar.SearchRangeX = [-InPar.SearchRangeX InPar.SearchRangeX];
end
if (numel(InPar.SearchRangeY)==1)
    InPar.SearchRangeY = [-InPar.SearchRangeY InPar.SearchRangeY];
end

% select best targets from List and Cat for matching
if (~isempty(InPar.FunCatSelect))
    FlagCat = InPar.FunCatSelect(Cat,InPar.FunCatSelectPar{:});
    Cat     = Cat(FlagCat,:);
end
if (~isempty(InPar.FunRefSelect))
    FlagRef = InPar.FunRefSelect(Ref,InPar.FunRefSelectPar{:});
    Ref     = Cat(FlagRef,:);
end

% % select best targets from List and Cat for matching
% % select targets if columns are in required range
% FlagSelCat = Util.array.array_select(Cat,InPar.CatSelectOp,InPar.CatSelect{:});
% FlagSelRef = Util.array.array_select(Ref,InPar.RefSelectOp,InPar.RefSelect{:});
% SubCat     = Cat(FlagSelCat,:);
% SubRef     = Ref(FlagSelRef,:);

    
% number of elements in Ref and Cat
Ncat = size(Cat,1);
Nref = size(Ref,1);

% if InPar.MinNinPeak<1 then use it as the fraction of matches for Ncat
if (InPar.MinNinPeak<1)
    InPar.MinNinPeak = ceil(InPar.MinNinPeak.*Ncat);
end


%Xref = Ref(:,InPar.ColXr);
%Yref = Ref(:,InPar.ColYr);

% if (InPar.MaxPeaksCheck>0 && ~issorted(Ycat))
%     error('Input catalog is not sorted');
% end

Nflip = size(InPar.Flip,1);
Res = [];
for Iflip=1:1:Nflip
    % Vectors of X/Y coordinates
    Xcat = Cat(:,InPar.ColXc);
    Ycat = Cat(:,InPar.ColYc);

    % Flip coordinates
    % Note this assumes that the coordinates are specified relative to
    % image center.
    Xcat = InPar.Flip(Iflip,1).*Xcat;
    Ycat = InPar.Flip(Iflip,2).*Ycat;
    

    Nrot = numel(Rotation);
    for Irot=1:1:Nrot
        % rotate Cat - rotation in deg
        % note the rotation is applied after flip!
        XcatRot = Xcat.*cosd(Rotation(Irot)) - Ycat.*sind(Rotation(Irot));
        YcatRot = Xcat.*sind(Rotation(Irot)) + Ycat.*cosd(Rotation(Irot));
        
        if (InPar.CutRefCat)
            MinX = min(XcatRot) - max(abs(InPar.SearchRangeX));
            MaxX = max(XcatRot) + max(abs(InPar.SearchRangeX));
            MinY = min(YcatRot) - max(abs(InPar.SearchRangeY));
            MaxY = max(YcatRot) + max(abs(InPar.SearchRangeY));

            % use only Ref sources in SubCat region
            FinCut = Ref(:,InPar.ColXr)>MinX & ...
                     Ref(:,InPar.ColXr)<MaxX & ...
                     Ref(:,InPar.ColYr)>MinY & ...
                     Ref(:,InPar.ColYr)<MaxY;
                 
            RefSub = Ref(FinCut,:);
        else
            RefSub = Ref;
        end

        [ResT,IndBest,H2]=ImUtil.pattern.match_pattern_shift([XcatRot,YcatRot],RefSub,...
                                                 'ColXc',1,...
                                                 'ColYc',2,...
                                                 'ColXr',InPar.ColXr,...
                                                 'ColYr',InPar.ColYr,...
                                                 'SearchRangeX',InPar.SearchRangeX,...
                                                 'SearchRangeY',InPar.SearchRangeY,...
                                                 'SearchRangeFactor',InPar.SearchRangeFactor,...
                                                 'SearchStepX',InPar.SearchStepX,...
                                                 'SearchStepY',InPar.SearchStepY,...
                                                 'MaxMethod',InPar.MaxMethod,...
                                                 'FalseAlarm',InPar.FalseAlarm,...
                                                 'RetSolFracMax',InPar.RetSolFracMax,...
                                                 'FunCatSelect',[],...
                                                 'FunCatSelectPar',{},...
                                                 'FunRefSelect',[],...
                                                 'FunRefSelectPar',{},...
                                                 'RetSolFracMax',InPar.RetSolFracMax,...
                                                 'MaxPeaks',InPar.MaxPeaks,...
                                                 'Radius',InPar.Radius,...
                                                 'SelectBest',InPar.SelectBest,...
                                                 'MaxPeaksCheck',0,...
                                                 'MinNinPeak',InPar.MinNinPeak);
        
        if (~isempty(IndBest))
            Npeak   = numel(ResT);
            CellRot = num2cell(Rotation(Irot).*ones(Npeak,1));
            [ResT(1:1:Npeak).Rot]   = deal(CellRot{:});
            CellFlipX = num2cell(InPar.Flip(Iflip,1).*ones(Npeak,1));
            CellFlipY = num2cell(InPar.Flip(Iflip,2).*ones(Npeak,1));
            [ResT(1:1:Npeak).FlipX]   = deal(CellFlipX{:});
            [ResT(1:1:Npeak).FlipY]   = deal(CellFlipY{:});
            if (isempty(Res))
                Res = ResT;
            else
                Res = [Res(:); ResT(:)];  
            end
        end
    end
end

% return only solutions that their HistNmatch> MaxHistNmatch.*InPar.RetSolFracMax
if (~isempty(Res))
    % sort by HistNmatch
    [~,SI] = sort([Res.HistNmatch],2,'descend');
    Res    = Res(SI);
    
    MaxHistNmatch = max([Res.HistNmatch]);
    FlagFracMax   = [Res.HistNmatch]>(MaxHistNmatch.*InPar.RetSolFracMax);
    Res           = Res(FlagFracMax);
    
    
end


% Check if Res contains the multiple candidates of the same solution
if (numel(Res)>1)
    DD = sqrt(([Res.HistShiftX] - [Res.HistShiftX]').^2 + ([Res.HistShiftY] - [Res.HistShiftY]').^2 );
    CritDist = sqrt(InPar.SearchStepX.^2 + InPar.SearchStepY.^2).*1.5;
    [DI,DJ] = find(triu(DD)>CritDist);
    % keep only one solution 
    Res = Res([1;DI]);
end

Npeak      = numel(Res);
NpeakCheck = min(InPar.MaxPeaksCheck,Npeak);   % max. number of peaks to check
%Res   = struct_def({'IndRef','IndCat','MatchedCat','MatchedRef','MatchedResid','StdResid','Std','MeanErr'},Npeak,1);
for Ipeak=1:1:NpeakCheck
    % given possible shift
    % matche Cat and Ref 
    
    Xcat = Cat(:,InPar.ColXc);
    Ycat = Cat(:,InPar.ColYc);
%     Xcat = Res(Ipeak).FlipX.*Xcat;
%     Ycat = Res(Ipeak).FlipY.*Ycat;

    %!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!
    %[ResMatch,OutCat]     = ImUtil.pattern.test_match(Cat,Ref,[Res(Ipeak).ShiftX, Res(Ipeak).ShiftY, Res(Ipeak).Rot./RAD],'Radius',InPar.Radius);
    %[ResMatch,OutCat,CatS]= ImUtil.pattern.test_match([Xcat, Ycat],Ref,...
    %                                                  [Res(Ipeak).ShiftX, Res(Ipeak).ShiftY, Res(Ipeak).Rot./RAD],...
    %                                                  'Radius',InPar.Radius,...
    %                                                  'Flip',[Res(Ipeak).FlipX, Res(Ipeak).FlipY]);
    
    %inserting the whole catalog instead of only X-Y columns
    [ResMatch,OutCat,CatS]= ImUtil.pattern.test_match(Cat,Ref,...
                                                      [Res(Ipeak).ShiftX, Res(Ipeak).ShiftY, Res(Ipeak).Rot./RAD],...
                                                      'ColXc',InPar.ColXc,...
                                                      'ColYc',InPar.ColYc,...
                                                      'Radius',InPar.Radius,...
                                                      'Flip',[Res(Ipeak).FlipX, Res(Ipeak).FlipY]);
    %!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!
                                                  
    Res(Ipeak).Nmatch     = numel(ResMatch.IndCat);
    %Res(Ipeak).MatchedCat = OutCat(ResMatch.IndCat,:);
    Res(Ipeak).MatchedCat = CatS(ResMatch.IndCat,:);
    %Res(Ipeak).MatchedCat = Cat(ResMatch.IndCat,1:2);
    Res(Ipeak).MatchedRef = Ref(ResMatch.IndRef,:);
    Res(Ipeak).ResidX     = ResMatch.ResidX;
    Res(Ipeak).ResidY     = ResMatch.ResidY;
    Res(Ipeak).Resid      = sqrt(ResMatch.ResidX.^2 + ResMatch.ResidY.^2);
    Res(Ipeak).StdX       = ResMatch.StdX;
    Res(Ipeak).StdY       = ResMatch.StdY;
    Res(Ipeak).Std        = std(Res(Ipeak).Resid);
    Res(Ipeak).Rstd       = Util.stat.rstd(Res(Ipeak).Resid);
    Res(Ipeak).MeanErr    = Res(Ipeak).Rstd./sqrt(Res(Ipeak).Nmatch);
end

  
if (Npeak==0)
    IndBest = [];
else
    if (NpeakCheck==0)
        % can use only MaxHistMatch
        [~,IndBest] = min([Res.HistNmatch]);
    else

        switch lower(InPar.SelectBest)
            case 'meanerr'
                [~,IndBest] = min([Res.MeanErr]);
            case 'std'
                [~,IndBest] = min([Res.Std]);
            case 'rstd'
                [~,IndBest] = min([Res.Rstd]);
            case 'n'
                [~,IndBest] = max([Res.Nmatch]);
            case 'comb'
                MaxNmatch = max([Res.Nmatch]);
                IndGood = find([Res.Nmatch]>(MaxNmatch.*0.8));
                [~,IndBest] = min([Res(IndGood).MeanErr]);
                IndBest = IndGood(IndBest);
            otherwise
                error('Unknown SelectBest option');
        end
    end
end



% % order Res by quality of solution
% switch lower(InPar.SortResBy)
%     case 'n'
%         [~,SI] = sort([Res.Nmatch]);
%     otherwise
%         error('Unknown SortResBy option');
% end
% Res = Res(SI);

% go over all possible peaks in the 2D histogram of X/Y distances
% and for each possibility search for matches in the complete
% catalogs
%Npeak      = numel(Res);
%NpeakCheck = min(InPar.MaxPeaksCheckInd,Npeak);   % max. number of peaks to check
%Res   = struct_def({'IndRef','IndCat','MatchedCat','MatchedRef','MatchedResid','StdResid','Std','MeanErr'},Npeak,1);

% % select best
% Res        = Res(1:NpeakCheck);
% for Ipeak=1:1:NpeakCheck
%     % given possible shift
%     % matche Cat and Ref 
%     
%     [ResMatch,OutCat]     = ImUtil.pattern.test_match(Cat,Ref,[Res(Ipeak).ShiftX, Res(Ipeak).ShiftY],'Radius',InPar.Radius);
%     Res(Ipeak).Nmatch     = numel(ResMatch.IndCat);
%     Res(Ipeak).MatchedCat = OutCat(ResMatch.IndCat,:);
%     Res(Ipeak).MatchedRef = Ref(ResMatch.IndRef,:);
%     Res(Ipeak).ResidX     = ResMatch.ResidX;
%     Res(Ipeak).ResidY     = ResMatch.ResidY;
%     Res(Ipeak).Resid      = sqrt(ResMatch.ResidX.^2 + ResMatch.ResidY.^2);
%     Res(Ipeak).StdX       = ResMatch.StdX;
%     Res(Ipeak).StdY       = ResMatch.StdY;
%     Res(Ipeak).Std        = std(Res(Ipeak).Resid);
%     Res(Ipeak).Rstd       = Util.stat.rstd(Res(Ipeak).Resid);
%     Res(Ipeak).MeanErr    = Res(Ipeak).Std./sqrt(Res(Ipeak).Nmatch);
% end
