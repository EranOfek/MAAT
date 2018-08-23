function [GoodInd,BadTimes,ExpTime]=find_badtimes(Input,varargin)
% Find bad times in X-ray event files by identyfing access of events.
% Package: AstroX
% Description: Given an X-ray event file, look for time ranges in which
%              the background is elevated above the mean background rate
%              by a given amount.
% Input  : - A matrix in which one of the columns is time of events, or
%            a vector of time of events.
%            Alternatievely, this can be an event file in FITS format.
%          * Arbitrary number of pairs of input arguments:...,key,val,...
%            The following keywords ara available:
%            'TimeCol'  - Index of time column in input matrix.
%                         Default is 1.
%            'MethodStd'- Method by which to estimate error on number of
%                         events in bin.
%                         Options are: 'sqrt'|'poisson'|'std'|'rstd'.
%                         Default is 'sqrt'.
%            'NinBin'   - Tyipcal number of events in bin to use in the
%                         identification of bad times. Default is 20.
%            'Nstd'     - Number of StD above the median count rate which
%                         to reject as high background. Default is 4.
% Output : - Vector of indices of good events (i.e., events which are not
%            during bad times).
%          - Matrix of bad time ranges, each row correspond to one range.
%            The first column for the start time and the second column for
%            the end time.
%          - Total exposure time of good times.
% Tested : Matlab 7.13
%     By : Eran O. Ofek                    Feb 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [GoodInd,BadTimes,ExpTime]=AstroX.find_badtimes('acisf00366N005_evt2.fits');
% Reliable: 2
%------------------------------------------------------------------------------


DefV.TimeCol             = 1;
DefV.MethodStd           = 'sqrt';
DefV.NinBin              = 100;
DefV.Nstd                = 4;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (ischar(Input))
    % load FITS binary table
    Cat = FITS.read_table(Input);
    
    TimeCol    = find(strcmpi(Cat.ColCell,'time'));
    EventTimes = Cat.Cat(:,TimeCol);

    % read Good Times from FITS
    TableGT     = fitsread(Input,'BinTable',2);
    GoodTimes   = [TableGT{1}, TableGT{2}];
   
else
    EventTimes = Input(:,InPar.TimeCol);
    GoodTimes  = GoodTimes;
    Col        = [];
    Table      = [];
end


% Total Integration
StartT     = min(EventTimes);
EndT       = max(EventTimes);
TotExpTime = EndT - StartT;
NtotEvt    = numel(EventTimes);

% select bin size such that the mean number of events per bin is
% InPar.NinBin
MeanEventRate = NtotEvt./TotExpTime;   % [events/s]
WindowSize    = InPar.NinBin./MeanEventRate;   % [s]
Nwin          = ceil(TotExpTime./WindowSize);
WindowSize    = TotExpTime./Nwin;
MeanEventPerWin = WindowSize.*MeanEventRate;

% bin
TimeEdges = (StartT:WindowSize:EndT)';
Nev       = histcounts(EventTimes,TimeEdges);

% Std estimate
switch lower(InPar.MethodStd)
    case 'std'
        Std = std(Nev);
    case 'rstd'
        Std = Util.stat.rstd(Nev');
    case 'sqrt'
        Std = sqrt(MeanEventPerWin);
    case 'poisson'
        Std = Util.stat.poissconf(MeanEventPerWin,'1');
        Std = Std(2) - MeanEventPerWin;
    otherwise
        error('Unknown MethodStd option');
end
        
FgoodWin   = Nev<(InPar.Nstd.*Std+MeanEventPerWin);
AllRanges  = [TimeEdges(1:end-1), TimeEdges(2:end)];
BadTimes   = AllRanges(~FgoodWin,:);


GoodInd = Util.array.find_ranges(EventTimes,AllRanges(FgoodWin,:));

ExpTime = TotExpTime - sum(BadTimes(:,2)-BadTimes(:,1));





