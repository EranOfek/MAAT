function [Err,ErrErr]=err_cl(Data,Conf)
% Numerical estimate of percentiles.
% Package: Util.stat
% Description: Given a vector of data points, calculate the lower and
%              upper bounds of an interval that contains a given
%              precentage (P) of the data. (1-P)/2 of the data points are
%              below and above the lower and upper bounds, respectively.
%              Remove NaNs prior to calculation.
% Input  : - Data vector.
%          - Lower and upper bounds of an interval that contains a given
%            precentage (P) of the data.
%            default is [0.6827; 0.9545; 0.9973]
%            error is calculate to each one of the levels.
% Output : - Matrix of confidence interval (CI), the first column for
%            left CI and the second column to right CI.
%            Each line for each confidence level.
%          - Estimate the standard errors in the confidence intervals.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Nov 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Err,ErrErr]=Util.stat.err_cl(randn(10000,1));
% Reliable: 1
%--------------------------------------------------------------------------
if (nargin==1)
   Conf = [0.6827;0.9545;0.9973];
elseif (nargin==2)
   % do nothing
else
   error('Illigal number of input arguments');
end

Data = Data(isnan(Data)==0);

N        = length(Data);
SortData = sort(Data);
%Mean     = mean(Data);
%Median   = median(Data);

BoundFrac = 0.5.*(1 - Conf);
NLower    = BoundFrac.*N;
NUpper    = (1-BoundFrac).*N;

if (N<2)
    Err = zeros(length(Conf),2).*NaN;
    ErrErr = zeros(length(Conf),2).*NaN;
else
   % take mean value between floor an ceil
   LowerConfPos = zeros(1,length(Conf));
   UpperConfPos = zeros(1,length(Conf));
   for I=1:1:length(Conf)
      if (floor(NLower(I))<1)
         LowerConfPos(I) = SortData(1);
      else
         LowerConfPos(I) = 0.5.*(SortData(floor(NLower(I))) + SortData(ceil(NLower(I))));
      end

      if (ceil(NUpper(I))>N)
         LowerConfPos(I) = SortData(N);
      else
         UpperConfPos(I) = 0.5.*(SortData(floor(NUpper(I))) + SortData(ceil(NUpper(I))));
      end
   end
   Err = [LowerConfPos', UpperConfPos'];
   
   if (nargout==2)
      % poisson noise
      N_Err = sqrt(NLower);
   
      LowerConfPosEU = zeros(1,length(Conf));
      LowerConfPosEL = zeros(1,length(Conf));
      for I=1:1:length(Conf)
         if ((ceil(NLower(I)+N_Err(I))>N) || (floor(NLower(I)+N_Err(I))<1))
            LowerConfPosEU = SortData(2);
         else
            LowerConfPosEU(I) = 0.5.*(SortData(floor(NLower(I)+N_Err(I))) + SortData(ceil(NLower(I)+N_Err(I))));
         end
      end
   
      for I=1:1:length(Conf)
         if (floor(NLower(I)-N_Err(I))<1)
            LowerConfPosEL(I) = SortData(1);
         else
            LowerConfPosEL(I) = 0.5.*(SortData(floor(NLower(I)-N_Err(I))) + SortData(ceil(NLower(I)-N_Err(I))));
         end
      end
      ErrErr = (LowerConfPosEU - LowerConfPosEL)';
   end
end
