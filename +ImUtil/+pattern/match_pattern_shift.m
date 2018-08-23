function [Res,IndBest,H2]=match_pattern_shift(Cat,Ref,varargin)
% Match X/Y coordinates in two lists in which the coordinates are shifted.
% Package: ImUtil.pattern
% Description: iven two lists containing two dimensional planar
%              coordinates (X, Y), find the possible shift transformations
%              between the two lists, using a histogram of all possible
%              combination of X and Y distances.         
% Input  : - Catalog list, containing at least two columns of [X,Y].
%            Should be sorted by Y.
%            If this is an AstCat object then AstCat/pattern_match_shift
%            will be called.
%          - Reference list, containing at least two columns of [X,Y].
%            Sorted by Y.
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
%                      'rstd'    - Choose the solution with the minmum
%                                  robust std.
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
% Example: [Res,IndBest,H2]=ImUtil.pattern.match_pattern_shift
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin==0)
    % simulation mode
    Nstar = 1000;
    Ref = rand(Nstar,2).*2048;
    Ref = sortrows(Ref,2);
    Noverlap = 100;
    Cat = [Ref(1:Noverlap,1), Ref(1:Noverlap,2)];
    Cat = [Cat; rand(Nstar-Noverlap,2).*2048];
    Cat(:,1) = Cat(:,1) + 20 + randn(Nstar,1).*0.3;
    Cat(:,2) = Cat(:,2) + 30 + randn(Nstar,1).*0.3;
    Cat      = sortrows(Cat,2);
end



DefV.ColXc            = 1;            % or string
DefV.ColYc            = 2;
DefV.ColXr            = 1;
DefV.ColYr            = 2;
DefV.FunCatSelect     = [];
DefV.FunCatSelectPar  = {};
DefV.FunRefSelect     = [];
DefV.FunRefSelectPar  = {};
DefV.SearchRangeX     = [-1000 1000];
DefV.SearchRangeY     = [-1000 1000];
DefV.SearchRangeFactor= 0.5;
DefV.SearchStepX      = 2;
DefV.SearchStepY      = 2;
DefV.MaxMethod        = 'adaptive'; %'sort';       % {'max1'|'sort'|'
DefV.FalseAlarm       = 1e-2;
DefV.MinNinPeak       = 5;            % min number of matches in peak - if <1 then fraction of ListCat
DefV.RetSolFracMax    = 0.5;          % return only solutions that their HistNmatch is larger than this fraction*max(HistNmatch)
DefV.MaxPeaks         = 10;           % maximum number of peaks to return
DefV.MaxPeaksCheck    = 10;           % out of the returned peaks - this is the maximum nuber of peaks to check
%DefV.CooType          = 'plane';
DefV.Radius           = 2;
%DefV.SearchMethod     = 'binms';
DefV.SelectBest       = 'comb';    % {'N','std','meanerr','comb'}

DefV.RegionMaxConn = 8;

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


if (InPar.MaxPeaksCheck>0 && ~issorted(Cat(:,InPar.ColYc)))
    Cat = sortrows(Cat,InPar.ColYc);
    %error('Input catalog is not sorted');
end

% Vectors of X/Y coordinates
Xcat = Cat(:,InPar.ColXc);
Ycat = Cat(:,InPar.ColYc);
Xref = Ref(:,InPar.ColXr);
Yref = Ref(:,InPar.ColYr);


% calculate matrices of X and Y distances
%Dx = bsxfun(@minus,Xref,Xcat.');
Dx = Xref - Xcat.';
%Dy = bsxfun(@minus,Yref,Ycat.');
Dy = Yref - Ycat.';

% generate a 2D histogram of X and Y distances
[H2,VecX,VecY] = Util.stat.hist2d(Dx(:),Dy(:),InPar.SearchRangeX,InPar.SearchStepX,InPar.SearchRangeY,InPar.SearchStepY);

switch lower(InPar.MaxMethod)
    case 'max1'
        InPar.MaxPeaks = 1;
end

StructFields = {'HistNmatch','HistShiftX','HistShiftY','HistPeakSN'};
Res = Util.struct.struct_def(StructFields,InPar.MaxPeaks,1);

% StD may sometime be zero due to rejection of all >0
% values in H2 - have to fix this
% meanwhile added max(std,1)...
ErrCL     = quantile(H2(:),0.99);
Std       = std(H2(H2<ErrCL));
Std       = max(1,Std);

% select peaks in 2D histogram of X and Y distances
% Res = struct('HistNmatch',[],'HistShiftX',[],'HistShiftY',[],'HistPeakSN',[],...
%              'ShiftX',[],'ShiftY',[],'Tran',[],...
%              'Rot',[],'FlipX',[],'FlipY',[]);
         
switch lower(InPar.MaxMethod)
    case 'max1'
        % select the highest maximum
        [MaxH,MaxInd] = max(H2(:));
        [MaxI,MaxJ]   = ind2sub(size(H2),MaxInd);
        
        %[MaxH,MaxIJ]    = Util.stat.maxnd(H2);
        
        if (MaxH>=InPar.MinNinPeak)
            Res.HistNmatch      = MaxH;
            Res.HistShiftX      = VecX(MaxJ);
            Res.HistShiftY      = VecY(MaxI);
            Res.HistPeakSN      = MaxH./Std;
            Nh                  = 1;
        else
            Res.HistNmatch      = NaN;
            Res.HistShiftX      = NaN;
            Res.HistShiftY      = NaN;
            Res.HistPeakSN      = NaN;
            Nh                  = 0;
        end
   
    case {'sort','adaptive'}
        switch lower(InPar.MaxMethod)
            case 'adaptive'
                % search peaks using matched filter
                
                DetThresh = poissinv(1-InPar.FalseAlarm,Std);
            case 'sort'
                DetThresh = InPar.MinNinPeak;
        end
        
%         tic;
%         H20 = H2;
%         H20(H20<DetThresh) = 0;
%         BW = imregionalmax(H20);
%         toc
         
        % Select N highest peaks
        % must use find as Ih2 is being used later on...
        Ih2 = find(H2(:)>=DetThresh);
        [SH,SI]         = sort(H2(Ih2),1,'descend');
        Nh              = numel(SH);
        % make sure that MaxPeaks is not larger then the number of elements
        InPar.MaxPeaks  = min(Nh,InPar.MaxPeaks);  
        SH              = SH(1:1:InPar.MaxPeaks);
        SI              = SI(1:1:InPar.MaxPeaks);
        
        % S/N of peak
        CellSN          = num2cell(SH(end-InPar.MaxPeaks+1:end)./Std);
        [Res(1:1:InPar.MaxPeaks).HistPeakSN]   = deal(CellSN{:});
        
        CellSH          = num2cell(SH(end-InPar.MaxPeaks+1:end));
        [Res(1:1:InPar.MaxPeaks).HistNmatch]    = deal(CellSH{:});
        
        [MaxI,MaxJ]     = ind2sub(size(H2),Ih2(SI));
        CellMaxJ        = num2cell(VecX(MaxJ));
        CellMaxI        = num2cell(VecY(MaxI));
        [Res(1:1:InPar.MaxPeaks).HistShiftX]   = deal(CellMaxJ{:});
        [Res(1:1:InPar.MaxPeaks).HistShiftY]   = deal(CellMaxI{:});
        
    otherwise
        error('Unknown MaxMethod option');
end

if (Nh==0)
    IndBest = [];
else

    % return only solutions that their HistNmatch> MaxHistNmatch.*InPar.RetSolFracMax
    MaxHistNmatch = max([Res.HistNmatch]);
    FlagFracMax   = [Res.HistNmatch]>(MaxHistNmatch.*InPar.RetSolFracMax);
    Res           = Res(FlagFracMax);

    Nt = numel(Res);

    % find exact centers of the peaks in the histogram
    for It=1:1:Nt
        Res(It).ShiftX = Res(It).HistShiftX;
        Res(It).ShiftY = Res(It).HistShiftY;
        Res(It).Tran   = affine2d([1 0 0; 0 1 0; Res(It).ShiftX, Res(It).ShiftY 1]);
    end


    % go over all possible peaks in the 2D histogram of X/Y distances
    % and for each possibility search for matches in the complete
    % catalogs
    Npeak      = numel(Res);
    NpeakCheck = min(InPar.MaxPeaksCheck,Npeak);   % max. number of peaks to check
    %Res   = struct_def({'IndRef','IndCat','MatchedCat','MatchedRef','MatchedResid','StdResid','Std','MeanErr'},Npeak,1);
    for Ipeak=1:1:NpeakCheck
        % given possible shift
        % matche Cat and Ref 

        [ResMatch,OutCat]     = ImUtil.pattern.test_match(Cat,Ref,[Res(Ipeak).ShiftX, Res(Ipeak).ShiftY],'Radius',InPar.Radius);
        Res(Ipeak).Nmatch     = numel(ResMatch.IndCat);
        Res(Ipeak).MatchedCat = OutCat(ResMatch.IndCat,:);
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
end