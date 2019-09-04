function PS=period_pdm(LC,FreqVec,varargin)
% SHORT DESCRIPTION HERE
% Package: timeseries
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: PS=timeseries.period_pdm([T,M,E],FreqVec);
% Reliable: 
%--------------------------------------------------------------------------


DefV.Nbin                 = 10;
DefV.ColT                 = 1;
DefV.ColM                 = 2;
DefV.ColE                 = 3;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% make sure FreqVec is a row vector
FreqVec = FreqVec(:).';

T = LC(:,InPar.ColT);
M = LC(:,InPar.ColM);
E = LC(:,InPar.ColE);

Nep   = numel(T);
Nfreq = numel(FreqVec);

% PhaseMatrix contains in each column the phase of the folded LC
% where each column corresponds to a different frequency
PhaseMatrix = mod(T,1./FreqVec).*FreqVec;

% replicate M
MatM = repmat(M,1,Nfreq);

% The index of the point in the histograms
IndHist = ceil(PhaseMatrix.*InPar.Nbin);

BinVar  = nan(InPar.Nbin,Nfreq);
BinMean = nan(InPar.Nbin,Nfreq);
BinMed  = nan(InPar.Nbin,Nfreq);
BinNpt  = nan(InPar.Nbin,Nfreq);

for Ibin=1:1:InPar.Nbin
    MatBins = nan(Nep,Nfreq);
    MatBins(IndHist==Ibin) = MatM(IndHist==Ibin);
    
    % calculate the properties in the current bin for all frequencies
    BinVar(Ibin,:)  = nanvar(MatBins,0,1);
    BinMean(Ibin,:) = nanmean(MatBins,1);
    BinMed(Ibin,:)  = nanmedian(MatBins,1);
    BinNpt(Ibin,:)  = sum(~isnan(MatBins),1);
    
end

VariancePerFrequency = sum(BinVar,1);

PS = [FreqVec.', VariancePerFrequency.'];

% take into account the number of points in each bin - the std probability
% distribution depands on the number of points

