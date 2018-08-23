function []=pattern_match_shift(Cat,Ref,varargin)
%--------------------------------------------------------------------------
% pattern_match_shift function                               class/@AstCat
% Description: Given two AstCat objects containing two dimensional planar
%              coordinates (X, Y), find the possible shift transformations
%              between the two lists, using a histogram of all possible
%              combination of X and Y distances. 
%              
% Input  : - An AstCat object. Each catalog contains at least two columns
%            of [X,Y]. Catalogs should be sorted by the Y column.
%            Each catalog will be matched agains the reference.
%          - An AstCat containing a single reference catalog. The catalog
%            contains at least two columns of [X,Y].
%            Catalog should be sorted by the Y column.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ColCatXY' - Cell array of column names or a vector of column
%                      indices of the X and Y coordinate in the catalog.
%                      Default is [1 2].
%            'ColRefXY' - Cell array of column names or a vector of column
%                      indices of the X and Y coordinate in the reference.
%                      Default is [1 2].
%            'CatSelect' - A string containing the arithmetic operation
%                      for selecting sources in the input catalog for
%                      matching. E.g., 'MAG_APER>14 & MAG_APER<18'.
%                       If empty, then select all entries.
%                       Default is empty.
%            'RefSelect' - A string containing the arithmetic operation
%                      for selecting sources in the input reference for
%                      matching. E.g., 'MAG_APER>14 & MAG_APER<18'.
%                       If empty, then select all entries.
%                       Default is empty.
%            'SearchRangeX' - X shift search range (in coordinate units).
%                      Default is [-1000 1000].
%            'SearchRangeY' - Y shift search range (in coordinate units).
%                      Default is [-1000 1000].
%            'SearchStepX' - X shift search step size (in coordinate units).
%                      Default is 1.
%            'SearchStepY' - Y shift search step size (in coordinate units).
%                      Default is 1.
%            'MaxMethod' - Method by which to use the best match
%                      solution/s. Options are:
%                      'sort' - sort by number of matches and return the
%                               best solutions.
%                      'max1' - Return the soltion with the largest
%                               number of matches.
%                      'adaptive' - Adapptive. Set the threshold using
%                               Posisson statistics and the required false
%                               alarm probability. Default.
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
%            'SearchRad' - Search radius for final matching. Default is 2.
%            'SearchMethod' - Catalog matching search method.
%                      See search_cat.m for options. Default is 'binms'.
%            'CheckSort' - Check if catalogs are sorted by Y {true|false}.
%                      Default is tru.
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




%            'ColXc' - Column name or column index of the X coordinate in
%                      the catalog. Default is 1.
%            'ColYc' - Column name or column index of the Y coordinate in
%                      the catalog. Default is 2.
%            'ColXr' - Column name or column index of the X coordinate in
%                      the reference. Default is 1.
%            'ColYr' - Column name or column index of the Y coordinate in
%                      the reference. Default is 2.
%            'CatSelect' - Cell array of column based selection criteria
%                      for catalog sources to match (e.g., {Col Min Max} or
%                      {Col @Fun} or {Col @Fun Val}).
%                      See array_select.m for details.
%                      Default is {}.
%            'CatSelectOp' - Catalog delection operator for array_select.m.
%                      Default is @all.
%            'RefSelect' - Cell array of column based selection criteria
%                      for reference sources to match (e.g., {Col Min Max} or
%                      {Col @Fun} or {Col @Fun Val}).
%                      See array_select.m for details.
%                      Default is {}.
%            'RefSelectOp' - Reference selection operator for array_select.m.
%                      Default is @all.
%            'SearchRangeX' - X shift search range (in coordinate units).
%                      Default is [-1000 1000].
%            'SearchRangeY' - Y shift search range (in coordinate units).
%                      Default is [-1000 1000].
%            'SearchStepX' - X shift search step size (in coordinate units).
%                      Default is 1.
%            'SearchStepY' - Y shift search step size (in coordinate units).
%                      Default is 1.
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
%          - The 2D histogram of ll possible combination of X and Y
%            distances.         
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------
RAD           = 180./pi;
%ARCSEC_IN_DEG = 3600;

CatField        = 'Cat';
ColCellField    = 'ColCell';
ColField        = 'Col';

DefV.ColCatXY      = [1 2];
DefV.ColRefXY      = [1 2];
DefV.CatSelect     = [];
DefV.RefSelect     = [];
% DefV.ColXc         = 1;            % or string
% DefV.ColYc         = 2;
% DefV.ColXr         = 1;
% DefV.ColYr         = 2;
% DefV.CatSelect     = {};             % {Col Min Max} or {Col @Fun} or {Col @Fun Val}
% DefV.CatSelectOp   = @all;
% DefV.RefSelect     = {};             % 
% DefV.RefSelectOp   = @all;
DefV.SearchRangeX  = [-1000 1000];
DefV.SearchRangeY  = [-1000 1000];
DefV.SearchStepX   = 2;
DefV.SearchStepY   = 2;
DefV.MaxMethod     = 'adaptive'; %'sort';       % {'max1'|'sort'|'
DefV.MinNinPeak    = 5;            % min number of matches in peak - if <1 then fraction of ListCat
DefV.MaxPeaks      = 10;           % maximum number of peaks to return
DefV.MaxPeaksCheck = 10;           % out of the returned peaks - this is the maximum nuber of peaks to check
%DefV.CooType       = 'plane';
DefV.SearchRad     = 2;
DefV.SearchMethod  = 'binms';
DefV.CheckSort     = true;
DefV.SelectBest    = 'comb';    % {'N','std','meanerr','comb'}
DefV.FalseAlarm    = 1e-4;
DefV.RegionMaxConn = 8;

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

% treat Ref
RefCol = colname2ind(Ref,InPar.ColRefXY);
% select Ref columns
RefFlag = col_arith(Ref, InPar.RefSelect);
Ref     = row_select(Ref,RefFlag);
% number of elements in Ref
Nref = size(Ref.(CatField),1);
Xref = Ref(:,RefCol(1));
Yref = Ref(:,RefCol(2));
% check sort
if (InPar.CheckSort),
    if (~issorted(Yref)),
        error('Reference catalog is not sorted');
    end
end



Ncat = numel(Cat);
% for each Cat
ResCat = struct('Res',cell(Ncat,1));
for Icat=1:1:Ncat,
    % Cat columns
    CatCol = colname2ind(Cat,InPar.ColCatXY);
    
    % select Cat columns
    CatFlag = col_arith(Cat(Icat), InPar.CatSelect);
    Cat(Icat)     = row_select(Cat(Icat),CatFlag);
    
    
    % number of elements in Cat
    Ncat = size(Cat(Icat).(CatField),1);

    % if InPar.MinNinPeak<1 then use it as the fraction of matches for Ncat
    if (InPar.MinNinPeak<1),
        InPar.MinNinPeak = ceil(InPar.MinNinPeak.*Ncat);
    end

    
    Xcat = Cat(:,CatCol(1));
    Ycat = Cat(:,CatCol(2));
    % check sort
    if (InPar.CheckSort),
        if (~issorted(Ycat)),
            error('Catalog is not sorted');
        end
    end
    
    % calculate matrices of X and Y distances
    Dx = bsxfun(@minus,Xref,Xcat.');
    Dy = bsxfun(@minus,Yref,Ycat.');
    % generate a 2D histogram of X and Y distances
    [H2,VecX,VecY] = hist2d(Dx(:),Dy(:),InPar.SearchRangeX,InPar.SearchStepX,InPar.SearchRangeY,InPar.SearchStepY);


    % select peaks in 2D histogram of X and Y distances
    switch lower(InPar.MaxMethod)
        case 'max1'
            % select the highest maximum
            [MaxH,MaxIJ]    = maxnd(H2);

            if (MaxH>=InPar.MinNinPeak),
                ResCat(Icat).Res.MaxHistMatch    = MaxH;
                ResCat(Icat).Res.MaxHistShiftX   = VecX(MaxIJ(2));
                ResCat(Icat).Res.MaxHistShiftY   = VecY(MaxIJ(1));
                Nh                  = 1;
            else
                ResCat(Icat).Res.MaxHistMatch    = NaN;
                ResCat(Icat).Res.MaxHistShiftX   = NaN;
                ResCat(Icat).Res.MaxHistShiftY   = NaN;
                Nh                  = 0;
            end
   
        case {'sort','adaptive'}
            switch lower(InPar.MaxMethod)
                case 'adaptive'
                    % search peaks using matched filter
                    ErrCL     = err_cl(H2(:),0.99);
                    Std       = std(H2(H2<ErrCL(2)));
                    DetThresh = poissinv(1-InPar.FalseAlarm,Std);
                case 'sort'
                    DetThresh = InPar.MinNinPeak;
            end
            % Select N highest peaks
            Ih2 = find(H2(:)>=DetThresh);
            [SH,SI]         = sort(H2(Ih2),1,'descend');

            Nh              = numel(SH);
            % make sure that MaxPeaks is not larger then the number of elements
            InPar.MaxPeaks  = min(Nh,InPar.MaxPeaks);  
            SH              = SH(1:1:InPar.MaxPeaks);
            SI              = SI(1:1:InPar.MaxPeaks);

            CellSH          = num2cell(SH(end-InPar.MaxPeaks+1:end));
            [ResCat(Icat).Res(1:1:InPar.MaxPeaks).MaxHistMatch]    = deal(CellSH{:});

            [MaxI,MaxJ]     = ind2sub(size(H2),Ih2(SI));
            CellMaxJ        = num2cell(VecX(MaxJ));
            CellMaxI        = num2cell(VecY(MaxI));
            [ResCat(Icat).Res(1:1:InPar.MaxPeaks).MaxHistShiftX]   = deal(CellMaxJ{:});
            [ResCat(Icat).Res(1:1:InPar.MaxPeaks).MaxHistShiftY]   = deal(CellMaxI{:});

        otherwise
            error('Unknown MaxMethod option');
    end
    if (Nh==0)
        ResCat(Icat).Res = [];
    end

    Nt = numel(ResCat(Icat).Res);

    % find exact centers of the peaks in the histogram
    for It=1:1:Nt,
        ResCat(Icat).Res(It).ShiftX = Res(It).MaxHistShiftX;
        ResCat(Icat).Res(It).ShiftY = Res(It).MaxHistShiftY;
        ResCat(Icat).Res(It).Tran   = affine2d([1 0 0; 0 1 0; ResCat(Icat).Res(It).ShiftX, ResCat(Icat).Res(It).ShiftY 1]);
    end

% transformPointsForward(TranRes.Tran,1,1)



% go over all possible peaks in the 2D histogram of X/Y distances
% and for each possibility search for matches in the complete
% catalogs
Npeak = min(InPar.MaxPeaksCheck,numel(ResCat(Icat).Res));   % max. number of peaks to check
%Res   = struct_def({'IndRef','IndCat','MatchedCat','MatchedRef','MatchedResid','StdResid','Std','MeanErr'},Npeak,1);
for Ipeak=1:1:Npeak,
    % given possible shift
    % matche Cat and Ref 
    SortedCat = Cat;
    % in the general case use e.g., transformPointsForward(TranRes.Tran,1,1)
    SortedCat(:,InPar.ColXc) = Cat(:,InPar.ColXc) + ResCat(Icat).Res(Ipeak).ShiftX;
    SortedCat(:,InPar.ColYc) = Cat(:,InPar.ColYc) + ResCat(Icat).Res(Ipeak).ShiftY;

    % This is not needed - assume cat are sorted!
    %SortedCat = sortrows(SortedCat,InPar.ColYc);

    ResMatch = search_cat(SortedCat(:,[InPar.ColXc, InPar.ColYc]),ResCat(Icat).Ref(:,[InPar.ColXr, InPar.ColYr]),[],...
                          'CooType',InPar.CooType,'SearchRad',InPar.SearchRad,'SearchMethod',InPar.SearchMethod);
                  
    ResCat(Icat).Res(Ipeak).Nmatch = sum([ResMatch.Nfound]==1);
    IndRef = find([ResMatch.Nfound]==1);
    IndCat = [ResMatch(IndRef).IndCat];
    ResCat(Icat).Res(Ipeak).IndRef       = IndRef;
    ResCat(Icat).Res(Ipeak).IndCat       = IndCat;
    ResCat(Icat).Res(Ipeak).MatchedCat   = SortedCat(IndCat,[InPar.ColXc, InPar.ColYc]);
    ResCat(Icat).Res(Ipeak).MatchedRef   = Ref(IndRef,[InPar.ColXr, InPar.ColYr]);
    ResCat(Icat).Res(Ipeak).MatchedResid = ResCat(Icat).Res(Ipeak).MatchedCat - ResCat(Icat).Res(Ipeak).MatchedRef;
    ResCat(Icat).Res(Ipeak).StdResid     = std(ResCat(Icat).Res(Ipeak).MatchedResid);
    if (ResCat(Icat).Res(Ipeak).Nmatch>1),
        ResCat(Icat).Res(Ipeak).Std          = sqrt(ResCat(Icat).Res(Ipeak).StdResid(1).^2 + ResCat(Icat).Res(Ipeak).StdResid(2).^2);

        % apply a grade for each solution
        ResCat(Icat).Res(Ipeak).MeanErr      = ResCat(Icat).Res(Ipeak).Std./sqrt(ResCat(Icat).Res(Ipeak).Nmatch);
    else
        % if only one match was found - std can't bes calculated
        ResCat(Icat).Res(Ipeak).Std = NaN;
        ResCat(Icat).Res(Ipeak).MeanErr = NaN;
    end
    % refine match using high order fit
    %Res(Ipeak).TranRes=fit_tran2d(Res(Ipeak).MatchedCat,Res(Ipeak).MatchedRef);
end


    if (Npeak==0),
        IndBest = [];
    else
        switch lower(InPar.SelectBest)
            case 'meanerr'
                [~,IndBest] = min([ResCat(Icat).Res.MeanErr]);
            case 'std'
                [~,IndBest] = min([ResCat(Icat).Res.Std]);
            case 'n'
                [~,IndBest] = max([ResCat(Icat).Res.Nmatch]);
            case 'comb'
                MaxNmatch = max([ResCat(Icat).Res.Nmatch]);
                IndGood = find([ResCat(Icat).Res.Nmatch]>(MaxNmatch.*0.8));
                [~,IndBest] = min([ResCat(Icat).Res(IndGood).MeanErr]);
                IndBest = IndGood(IndBest);
            otherwise
                error('Unknown SelectBest option');
        end
    end
end