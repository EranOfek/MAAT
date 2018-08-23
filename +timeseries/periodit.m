function [SuperP,SuperF,MidTime,FreqSer]=periodit(TimeSerMat,EpochPars,TaperPars,FreqVec,PeriodType,Disp)
% Calculate the periodogram as a function of time (will be removed)
% Package: timeseries
% Description: Calculate periodogram as function of time. This script can
%              be used for searchnig for long term varibility in
%              periodicities. The time series ins divided into blocks
%              with or without overlays. Each block can be tapered, and a 
%              periodogram (Clasical or Scargle) is calculated. The results
%              are displaied as a mesh of frequncy vs. time vs. power.
% Input  : - Matrix of observations. Should contain at least two
%            column. First column for time and second column for
%            "magnitude". Each row for each observation.
%          - vector of parameters defined the epochs of sub
%            periodograms. First element should contain the number
%            of epochs in the time span. and the second element
%            should contain the size of each sub epoch in units of
%            the total time span. default is [10,0.2].
%          - Tapering function parameters. Each sub epoch could be
%            weighted by a tapering function.
%            Avaliable tapers:
%            ['C',p] : cosine bell taper while "p" is the precentage
%                       of flat (uniform) part of the cosine bell in
%                       the range 0-1.
%            ['F']   : Flat (uniform) taper. No tapering.
%            default is ['C',0.6].
%          - frequency vector of the type:
%            [Low_Freqency, High_Frequency, Frequency_Interval].
%            default is [0,1./(2min(Interval)),1./(2(Time_Span)].
%          - Periodogram type:
%            'c' : classical periodogram.
%            's' : Scargle periodogram.
%            defualt is 'c' (classical).
%          - Display mesh at end of calculation.
%            'y' or 'n'. default is 'y'.
% Output : - Matrix of powers each column contain power for each
%            sub epoch.
%          - Matrix of leading frequency. Not Yet Avaliable.
%          - vector of sub epochs mid time.
%          - vector of frequencies searched.
% Remarks: Not all the option were checked. Please report
%          any bug/anomalies.
% Example: [P,F,T,Fr]=timeseries.periodit([(1:1:1000)',rand(1000,1)],[30,0.1],[],[],'s');
% Tested : Matlab 5.1
%     By : Eran O. Ofek                    Dec 1997
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------

ColT = 1;
ColO = 2;
[Nobs,Ncol]=size(TimeSerMat);
if Ncol<2,
   error('Number of columns in the TimeSerMat should be at least 2');
end

% check number of arguments.
if nargin==1,
   EpochPars  = [];
   TaperPars  = [];
   FreqVec    = [];
   PeriodType = 'c';   % periodogram default is classical (periodia.m).
   Disp       = 'y';   % default. display mesh.
elseif nargin==2,
   TaperPars  = [];
   FreqVec    = [];
   PeriodType = 'c';   % periodogram default is classical (periodia.m).
   Disp       = 'y';   % default. display mesh.
elseif nargin==3,
   FreqVec    = [];
   PeriodType = 'c';   % periodogram default is classical (periodia.m).
   Disp       = 'y';   % default. display mesh.
elseif nargin==4,
   PeriodType = 'c';   % periodogram default is classical (periodia.m).
   Disp       = 'y';   % default. display mesh.
elseif nargin==5,
   Disp       = 'y';   % default. display mesh.
else
   error('Number of input arguments should be 1,2,3,4 or 5');
end

TimeSpan = TimeSerMat(Nobs,ColT) - TimeSerMat(1,ColT); 

% handling the Epoch Parameters vector (EpochPars).
[EpochParsSizeI,EpochParsSizeJ]=size(EpochPars);
if (EpochParsSizeI==0),
   % null EpochPars
   % use defaults
   EpochNum    = 10;    % number of epochs in time span.
   EpochLength = 0.2;   % epoch length in units of time span.
else
   if (length(EpochPars)~=2),
      error('Number of elements in EpochPars vector should be 2');
   end
   EpochNum    = EpochPars(1);
   EpochLength = EpochPars(2);
end

% handling the tapering function parameters (TaperPars).
[TaperParsSizeI,TaperParsSizeJ]=size(TaperPars);
if (TaperParsSizeI==0),
   % null TaperPars
   % use defaults
   TaperFun  = 'C';   % cosine bell tapering function.
   TaperFunP = 0.6;    % tapering function parameter.
                       % in the case of 'C' this is the precentage
                       % of the flat part.
else
   if (length(TaperPars)==1),
      if TaperPars(1)=='F',
         TaperFun = 'F';     % Flat (uniform) tapering function;
      else
         error('Unknown Tapering function');
      end
      % add others single parameters tapering function
   elseif (length(TaperPars)==2),
      if TaperPars(1)=='C',
         TaperFun = 'C';     % cosine bell tapering function;
      else
         error('Unknown Tapering function');
      end
      % add others single parameters tapering function
   else
      error('Number of elements in TaperPars vector should be null, 1 or 2');
   end
end

% handling the Frequncy Vector (FreqVec)
[FreqVecSizeI,FreqVEcSizeJ]= size(FreqVec);
if (FreqVecSizeI==0),
   % null FreqVec
   % use defaults
   LowFreq  = 0;
   HighFreq = 1./(2.*min(abs(diff(TimeSerMat(:,ColT)))));
   FreqInt  = 1./(2.*TimeSpan);
else
   if (length(FreqVec)~=3),
      error('Number of elements in FreqVec should be null or 3');
   else
      LowFreq  = FreqVec(1);
      HighFreq = FreqVec(2);
      FreqInt  = FreqVec(3);
   end
end


% start program
StartEpoch = TimeSerMat(1,ColT);
EpochLenT  = EpochLength.*TimeSpan;
EpochStep  = TimeSpan./EpochNum;

%check P size

if (PeriodType=='c'),
   [P,F] = timeseries.periodia(TimeSerMat(1:ceil(0.1.*Nobs),:),LowFreq,HighFreq,FreqInt);
elseif (PeriodType=='s'),
   [P,F] = timeseries.periodis(TimeSerMat(1:ceil(0.1.*Nobs),:),LowFreq,HighFreq,FreqInt);
else
   error('unknown periodogram type');
end
FreqSer    = P(:,1);
SuperP     = zeros(length(P(:,1)),EpochNum-1);
SuperF     = zeros(10,EpochNum-1);
MidTime    = zeros(1,EpochNum-1);
for Ind=0:1:EpochNum-2,
   JJ = find(TimeSerMat(:,ColT) >= (StartEpoch+Ind.*EpochStep) & TimeSerMat(:,ColT) < (StartEpoch+Ind.*EpochStep+EpochLenT));
   NewMat = TimeSerMat(JJ,:);
   %NewMat = rangebrk(TimeSerMat, StartEpoch+Ind.*EpochStep, StartEpoch+Ind.*EpochStep+EpochLenT,ColT);
   if (length(NewMat(1,:))==1),
      % NewMat is empty. continue.
      MidTime(1,Ind+1) = StartEpoch+Ind.*EpochStep + 0.5.*EpochLenT;
   else
      MidTime(1,Ind+1) = StartEpoch+Ind.*EpochStep + 0.5.*EpochLenT;
      % tapering
      [WeightFun,Range] = timeseries.cosbell(TaperFunP,NewMat(:,ColT));
      NewMat(:,ColO) = NewMat(:,ColO).*WeightFun;
      if (PeriodType=='c'),
         [P,F] = timeseries.periodia(NewMat,LowFreq,HighFreq,FreqInt);
      elseif (PeriodType=='s'),
         [P,F] = timeseries.periodis(NewMat,LowFreq,HighFreq,FreqInt);
      else
         error('unknown periodogram type');
      end
      SuperP(:,Ind+1) = P(:,2);
      %SuperF(:,Ind+1) = F(length(F(:,1))-9:length(F(:,1)),1);
   end
end

if (Disp=='y'),
   mesh(MidTime,FreqSer,SuperP);
   xlabel('Time')
   ylabel('Frequency')
   zlabel('Power')
end

