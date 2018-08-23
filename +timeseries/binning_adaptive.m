function [B,OutCol]=binning_adaptive(Data,StartEnd,OutCol,varargin)
% Adaptive binning of time series. Bin size according to std in bin.
% Package: timeseries
% Description: Adaptive binning of a time series, where the bin size is
%              determined by the std of data in bin.
% Input  : - Data - e.g., [Time, Mag]
%          - [Start End] time. If empty, then use min and max of time.
%            Default is [].
%          - Cell array of column names to calculate. Possible column
%            names are:
%            'MidBin'
%            'StartBin'
%            'EndBin'
%            'MeanBin'
%            'MedianBin'
%            'StdBin'
%            Or any function that returns a scalar to apply to the
%            observations (e.g., @numel, @mean, @median).
%            There are several special functions that will be applied also
%            for the errors column. These are @wmean, @wmedian.
%            Default is {'MidBin', @numel, @mean, @median, @std}.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ColT' - Time column. Default is 1.
%            'ColM' - Mag column. Default is 2.
%            'TrialBin' - Vector of bin size to try.
%                     Default is [1 30 90].
%            'StdThresh' - Maximum std above to switch to latrger bin size.
%                     Default is 0.1.
%            'Frac' - If fraction of bins with either 0|1 points or std
%                     larger than StdThresh is larger than this fraction,
%                     then switch to larger bin. Default is 0.5.
%            'GapDef' - Time gap size above to sepearate the data.
%                     Default is 60.
%            'Clean' - Clean NaN values. Default is true.
% Output : - Matrix or an AstCat object with with the requested columns.
%          - Output column names.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [B,OutCol]=timeseries.binning_adaptive(Data,[],{'MidBin', @numel, @mean, @median, @Util.stat.rstd});
% Reliable: 
%--------------------------------------------------------------------------


Def.StartEnd = [];
Def.OutCol   = {'MidBin', @numel, @mean, @median, @std};
if (nargin==1)
    StartEnd = Def.StartEnd;
    OutCol   = Def.OutCol;
elseif (nargin==2)
    OutCol   = Def.OutCol;
else
    % do nothing
end
    

DefV.ColT                 = 1;
DefV.ColM                 = 2;
DefV.TrialBin             = [2 15 30 90];
DefV.StdThresh            = 0.05;
DefV.Frac                 = 0.5;
DefV.GapDef               = 60;
DefV.Clean                = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (isempty(StartEnd))
    StartEnd = [min(Data(:,InPar.ColT)), max(Data(:,InPar.ColT))];
end
    
Data = sortrows(Data,InPar.ColT);


Ntrial = numel(InPar.TrialBin);
Itrial = 1;

Time = Data(:,InPar.ColT);
DeltaT = diff([-Inf; Time; Inf]);

IndGap = find(DeltaT>InPar.GapDef);
Ngap = numel(IndGap);
for Igap=1:1:Ngap-1
    I1 = IndGap(Igap);
    I2 = IndGap(Igap+1)-1;
    
    Cont = true;
    Itrial = 1;
    while (Cont && Itrial<=Ntrial)
        BinT(Igap).BT = timeseries.binning(Data(I1:I2,:),InPar.TrialBin(Itrial),StartEnd,{@numel,@median,@Util.stat.rstd});
        BinT(Igap).B  = timeseries.binning(Data(I1:I2,:),InPar.TrialBin(Itrial),StartEnd,OutCol);
        BinT(Igap).FlagGood = BinT(Igap).BT(:,1)~=0 & BinT(Igap).BT(:,1)>1 & abs(BinT(Igap).BT(:,3)./BinT(Igap).BT(:,2))<InPar.StdThresh;
        BinT(Igap).FlagAll  = BinT(Igap).BT(:,1)~=0;

        Frac = sum(BinT(Igap).FlagGood)./sum(BinT(Igap).FlagAll);

        if (Frac<InPar.Frac)
            % fraction of good points is too small
            % switch to lager bins
            Cont = true;
        else
            Cont = false;
        end
        Itrial = Itrial + 1;
    end
    
    BinT(Igap).B = BinT(Igap).B.';
end    

B = [BinT.B].';
% clean B
if (InPar.Clean)
    B = B(~any(isnan(B),2),:);
end





