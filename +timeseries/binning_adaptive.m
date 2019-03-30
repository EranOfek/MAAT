function AllEdges=binning_adaptive(T,varargin)
% Construct list of edges for binning, given some bin size criteria
% Package: timeseries
% Description: 
% Input  : - List of times for which to construct edges for binning.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AllEdges=timeseries.binning_adaptive
% Reliable: 
%--------------------------------------------------------------------------


%T=[1 1.1   2 2.2 2.3 2.6 2.8 3 3.4 3.8 4 4.4  6 6.1]


DefV.MaxBin               = 1;
DefV.MinGap               = 0.5;
DefV.BreakMethod          = 'num'; % 'num' | 'time'
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

T  = T(:);
T  = sort(T);
DT = [diff(T);0];


IG = [0; find(DT>InPar.MinGap); numel(T)];
Ng = numel(IG)-1;
AllEdges = [];
for Ig=1:1:Ng
    Ind = ((IG(Ig)+1):1:IG(Ig+1)).';
    Range = range(T(Ind));
    NsubBins = ceil(Range./InPar.MaxBin);
    
    if (NsubBins>1)
        % break into several bins
        switch lower(InPar.BreakMethod)
            case 'num'
                NsubObs   = numel(Ind);
                NobsInBin = ceil(NsubObs./NsubBins);
                Edges     = T(Ind(1:NobsInBin:NsubObs));
            case 'time'
                SmallBinSize = Range./NsubBins;
                Edges = [min(T(Ind)):SmallBinSize:max(T(Ind))].';
                Edges = Edges(1:end-1);
            otherwise
                error('Unknown BreakMethod option');
        end
    else
        Edges = min(T(Ind));
    end
    AllEdges = [AllEdges;Edges];
end
AllEdges = [AllEdges;max(T)];
    


