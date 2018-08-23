function [Res,IndBest,H2]=pattern_match_shift(Cat,Ref,varargin)
% Match two [X,Y] catalogs by X/Y shift and pattern recognition 
% Package: ImUtil.Im
% Description: Given two lists containing two dimensional planar
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
%            'CatSelect' - Cell array of column based selection criteria
%                      for catalog sources to match (e.g., {Col Min Max} or
%                      {Col @Fun} or {Col @Fun Val}).
%                      See Util.array.array_select.m for details.
%                      Default is {}.
%            'CatSelectOp' - Catalog delection operator for
%                      Util.array.array_select.m.
%                      Default is @all.
%            'RefSelect' - Cell array of column based selection criteria
%                      for reference sources to match (e.g., {Col Min Max} or
%                      {Col @Fun} or {Col @Fun Val}).
%                      See Util.array.array_select.m for details.
%                      Default is {}.
%            'RefSelectOp' - Reference selection operator for
%                      Util.array.array_select.m.
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
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Out = FITS.read_table('*p100019*.ctlg');
%          [Res,IndBest]=ImUtil.Im.pattern_match_shift(Out(1).Cat(:,[3 4]),Out(2).Cat(:,[3 4]));
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin==0)
    % simulation mode
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



DefV.ColXc         = 1;            % or string
DefV.ColYc         = 2;
DefV.ColXr         = 1;
DefV.ColYr         = 2;
DefV.CatSelect     = {};             % {Col Min Max} or {Col @Fun} or {Col @Fun Val}
DefV.CatSelectOp   = @all;
DefV.RefSelect     = {};             % 
DefV.RefSelectOp   = @all;
DefV.SearchRangeX  = [-1000 1000];
DefV.SearchRangeY  = [-1000 1000];
DefV.SearchStepX   = 2;
DefV.SearchStepY   = 2;
DefV.MaxMethod     = 'adaptive'; %'sort';       % {'max1'|'sort'|'
DefV.MinNinPeak    = 5;            % min number of matches in peak - if <1 then fraction of ListCat
DefV.MaxPeaks      = 10;           % maximum number of peaks to return
DefV.MaxPeaksCheck = 10;           % out of the returned peaks - this is the maximum nuber of peaks to check
DefV.CooType       = 'plane';
DefV.SearchRad     = 2;
DefV.SearchMethod  = 'binms';
DefV.SelectBest    = 'comb';    % {'N','std','meanerr','comb'}
DefV.FalseAlarm    = 1e-4;
DefV.RegionMaxConn = 8;

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% deal with struct/SIM input
%if (isstruct(Cat) || SIM.issim(Cat)),    
%    [Cat,InPar.ColXc,InPar.ColYc]=simcat2array(Cat,{InPar.ColXc, InPar.ColYc});
%end
%if (isstruct(Ref) || SIM.issim(Ref)),
%    [Ref,InPar.ColXr,InPar.ColYr]=simcat2array(Ref,{InPar.ColXr, InPar.ColYr});
%end

    

% If SearchRange is a scalar then assume the range is symmetric
if (numel(InPar.SearchRangeX)==1)
    InPar.SearchRangeX = [-InPar.SearchRangeX InPar.SearchRangeX];
end
if (numel(InPar.SearchRangeY)==1)
    InPar.SearchRangeY = [-InPar.SearchRangeY InPar.SearchRangeY];
end


% select best targets from List and Cat for matching
% select targets if columns are in required range
FlagSelCat = Util.array.array_select(Cat,InPar.CatSelectOp,InPar.CatSelect{:});
FlagSelRef = Util.array.array_select(Ref,InPar.RefSelectOp,InPar.RefSelect{:});
SubCat     = Cat(FlagSelCat,:);
SubRef     = Ref(FlagSelRef,:);

    
% number of elements in Ref and Cat
Ncat = size(SubCat,1);
Nref = size(SubRef,1);

% if InPar.MinNinPeak<1 then use it as the fraction of matches for Ncat
if (InPar.MinNinPeak<1)
    InPar.MinNinPeak = ceil(InPar.MinNinPeak.*Ncat);
end

% Vectors of X/Y coordinates
Xcat = SubCat(:,InPar.ColXc);
Ycat = SubCat(:,InPar.ColYc);
Xref = SubRef(:,InPar.ColXr);
Yref = SubRef(:,InPar.ColYr);

if (InPar.MaxPeaksCheck>0 && ~issorted(Ycat))
    error('Input catalog is not sorted');
end

% calculate matrices of X and Y distances
Dx = bsxfun(@minus,Xref,Xcat.');
Dy = bsxfun(@minus,Yref,Ycat.');
% generate a 2D histogram of X and Y distances
[H2,VecX,VecY] = Util.stat.hist2d(Dx(:),Dy(:),InPar.SearchRangeX,InPar.SearchStepX,InPar.SearchRangeY,InPar.SearchStepY);

% select peaks in 2D histogram of X and Y distances
switch lower(InPar.MaxMethod)
    case 'max1'
        % select the highest maximum
        [MaxH,MaxIJ]    = maxnd(H2);
        
        if (MaxH>=InPar.MinNinPeak)
            Res.MaxHistMatch    = MaxH;
            Res.MaxHistShiftX   = VecX(MaxIJ(2));
            Res.MaxHistShiftY   = VecY(MaxIJ(1));
            Nh                  = 1;
        else
            Res.MaxHistMatch    = NaN;
            Res.MaxHistShiftX   = NaN;
            Res.MaxHistShiftY   = NaN;
            Nh                  = 0;
        end
   
    case {'sort','adaptive'}
        switch lower(InPar.MaxMethod)
            case 'adaptive'
                % search peaks using matched filter
                %ErrCL     = err_cl(H2(:),0.99);
                %Std       = std(H2(H2<ErrCL(2)));
                
                % StD may sometime be zero due to rejection of all >0
                % values in H2 - have to fix this
                % meanwhile added max(std,1)...
                ErrCL     = quantile(H2(:),0.99);
                Std       = std(H2(H2<ErrCL));
                Std       = max(1,Std);
                DetThresh = poissinv(1-InPar.FalseAlarm,Std);
            case 'sort'
                DetThresh = InPar.MinNinPeak;
        end
        % Select N highest peaks
        % must use find as Ih2 is being used later on...
        Ih2 = find(H2(:)>=DetThresh);
        [SH,SI]         = sort(H2(Ih2),1,'descend');
        Nh              = numel(SH);
        % make sure that MaxPeaks is not larger then the number of elements
        InPar.MaxPeaks  = min(Nh,InPar.MaxPeaks);  
        SH              = SH(1:1:InPar.MaxPeaks);
        SI              = SI(1:1:InPar.MaxPeaks);
        
        CellSH          = num2cell(SH(end-InPar.MaxPeaks+1:end));
        [Res(1:1:InPar.MaxPeaks).MaxHistMatch]    = deal(CellSH{:});
        
        [MaxI,MaxJ]     = ind2sub(size(H2),Ih2(SI));
        CellMaxJ        = num2cell(VecX(MaxJ));
        CellMaxI        = num2cell(VecY(MaxI));
        [Res(1:1:InPar.MaxPeaks).MaxHistShiftX]   = deal(CellMaxJ{:});
        [Res(1:1:InPar.MaxPeaks).MaxHistShiftY]   = deal(CellMaxI{:});
        
    otherwise
        error('Unknown MaxMethod option');
end
if (Nh==0)
    Res = [];
end

Nt = numel(Res);

% find exact centers of the peaks in the histogram
for It=1:1:Nt
    Res(It).ShiftX = Res(It).MaxHistShiftX;
    Res(It).ShiftY = Res(It).MaxHistShiftY;
    Res(It).Tran   = affine2d([1 0 0; 0 1 0; Res(It).ShiftX, Res(It).ShiftY 1]);
end

% transformPointsForward(TranRes.Tran,1,1)



% go over all possible peaks in the 2D histogram of X/Y distances
% and for each possibility search for matches in the complete
% catalogs
Npeak    = numel(Res);
NpeakCheck = min(InPar.MaxPeaksCheck,Npeak);   % max. number of peaks to check
%Res   = struct_def({'IndRef','IndCat','MatchedCat','MatchedRef','MatchedResid','StdResid','Std','MeanErr'},Npeak,1);
for Ipeak=1:1:NpeakCheck
    % given possible shift
    % matche Cat and Ref 
    SortedCat = Cat;
    % in the general case use e.g., transformPointsForward(TranRes.Tran,1,1)
    SortedCat(:,InPar.ColXc) = Cat(:,InPar.ColXc) + Res(Ipeak).ShiftX;
    SortedCat(:,InPar.ColYc) = Cat(:,InPar.ColYc) + Res(Ipeak).ShiftY;

    % This is not needed - assume cat are sorted!
    %SortedCat = sortrows(SortedCat,InPar.ColYc);

    ResMatch = VO.search.search_cat(SortedCat(:,[InPar.ColXc, InPar.ColYc]),Ref(:,[InPar.ColXr, InPar.ColYr]),[],...
                          'CooType',InPar.CooType,'SearchRad',InPar.SearchRad,'SearchMethod',InPar.SearchMethod);
                  
    Res(Ipeak).Nmatch = sum([ResMatch.Nfound]==1);
    IndRef = ([ResMatch.Nfound]==1);
    IndCat = [ResMatch(IndRef).IndCat];
    Res(Ipeak).IndRef       = IndRef;
    Res(Ipeak).IndCat       = IndCat;
    Res(Ipeak).MatchedCat   = SortedCat(IndCat,[InPar.ColXc, InPar.ColYc]);
    Res(Ipeak).MatchedRef   = Ref(IndRef,[InPar.ColXr, InPar.ColYr]);
    Res(Ipeak).MatchedResid = Res(Ipeak).MatchedCat - Res(Ipeak).MatchedRef;
    Res(Ipeak).StdResid     = std(Res(Ipeak).MatchedResid);
    if (Res(Ipeak).Nmatch>1)
        Res(Ipeak).Std          = sqrt(Res(Ipeak).StdResid(1).^2 + Res(Ipeak).StdResid(2).^2);

        % apply a grade for each solution
        Res(Ipeak).MeanErr      = Res(Ipeak).Std./sqrt(Res(Ipeak).Nmatch);
    else
        % if only one match was found - std can't bes calculated
        Res(Ipeak).Std = NaN;
        Res(Ipeak).MeanErr = NaN;
    end
    % refine match using high order fit
    %Res(Ipeak).TranRes=fit_tran2d(Res(Ipeak).MatchedCat,Res(Ipeak).MatchedRef);
end


if (Npeak==0)
    IndBest = [];
else
    if (NpeakCheck==0)
        % can use only MaxHistMatch
        [~,IndBest] = min([Res.MaxHistMatch]);
    else
        
        switch lower(InPar.SelectBest)
            case 'meanerr'
                [~,IndBest] = min([Res.MeanErr]);
            case 'std'
                [~,IndBest] = min([Res.Std]);
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