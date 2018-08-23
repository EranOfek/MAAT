function PDM=pdm(Data,FreqVec,Nbin)
% Periodicity search using period dispersion minimization
% Package: timeseries
% Description: Periodicity search using phase dispersion minimization.
%              The time series is folded into trail frequencies, and in
%              each trial frequency, the dispersion in bins is calculated.
% Input  : - Data matrix, first column for time and second for magnitude.
%          - Frequency vector, [minimum_freq, maximum_freq, delta_freq]
%	   - Number of bins in the folded period.
% Output : - Matrix in which the first column is frequency,
%            the second is the dispersion indicator normalized
%            by the variance of the data (for a true period =1).
% Reference: Stellingwerf, R.F. ApJ 224, 953-960 (1978).
%            Schwarzenberg-Czerny, A. ApJL, 489, 941-945 (1997)
% see also : minds function.
% Tested : Matlab 5.0
%     By : Eran O. Ofek                    Jun 1994
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: t=(1:1:1000)'; M=sin(2.*pi.*t./27.1)+randn(1000,1).*0.1; E=0.1.*ones(1000,1); P=timeseries.pdm([t, M E],(0:0.0005:1)',10);
% Reliable: 
%------------------------------------------------------------------------------

Col.T   = 1;
Col.M   = 2;
Col.E   = 3;

Def.Err = 1;

Ndata = size(Data,1);
if (size(Data,2)==2)
    Data = [Data, ones(Ndata,1).*Def.Err];
end


% null hypothesis
Mean    = mean(Data(:,Col.M));
Std     = std(Data(:,Col.M));
Chi2_H0 = sum(((Data(:,Col.M) - Mean)./Data(:,Col.E)).^2);

BinSize = 1./Nbin;

%Var = std(Data(:,2),0).^2;
%Ntot = length(Data(:,1));
%D0 = Ntot - 1;
%D1 = BinNum - 1;
%D2 = Ntot - BinNum;

%PDM = zeros(size([FreqVec(1):FreqVec(3):FreqVec(2)],2));

%I = 0;

Nf = numel(FreqVec);
PDM.Freq = FreqVec;
PDM.DeltaChi2 = zeros(Nf,1);
PDM.Chi2 = zeros(Nf,1);
PDM.Dof  = zeros(Nf,1);
PDM.rms  = zeros(Nf,1);
for If=1:1:Nf
    
   Phase = [Data(:,Col.T).*FreqVec(If) - fix(Data(:,Col.T).*FreqVec(If))];
   BinInd = ceil(Phase./BinSize);
   BinInd(BinInd==0) = 1;
   
   B      = timeseries.binning([Phase,Data(:,Col.M)],BinSize,[0 1],{'MeanBin',@mean,@std,@numel});
   
   Chi2   = sum(((Data(:,Col.M) - B(BinInd,2))./Data(:,Col.E)).^2);
   
   PDM.DeltaChi2(If) = Chi2 - Chi2_H0;
   PDM.Chi2(If)      = Chi2;
   PDM.Dof(If)       = Nbin;
   PDM.rms(If)       = std(Data(:,Col.M) - B(BinInd,2));
end

   
%    
%    
%    [BFData] = timeseries.binning_old(Folded,1./BinNum,0,1);
% 
%    N   = BFData(:,7);
%    StD = BFData(:,3).*sqrt(BFData(:,7));
%    S2  = sum((N-1).*StD.^2)./(sum(N) - BinNum);
%    Theta = S2./Var;
%    PDM(I,1:2) = [Freq, Theta]
