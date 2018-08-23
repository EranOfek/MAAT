function S=perioent(x,FreqVec,m)
% Periodogram using information entropy
% Package: timeseries
% Description: Periodogram using information entropy. For each trail
%              period, the phase-magnitude space is divided into m by m
%              cells, and the information entropy is calculated.
%              Then the probability to find observations in each square is
%              MU(i) and the Entropy is (-Sigma{MU*ln(MU)}). The output
%              Entropy is normalized by ln(m*m).
% Input  : - Data matrix, sorted by time.
%          - Frequency Vector, [Lower Freq., Higher Freq., Freq. Interval]
%          - Square root of number of elements in the unit square. 
% Output : - Two column matrix, [Frequency, Entropy].
% reference: Cincotta, Mendez & Nunez 1995, ApJ 449, 231-235.
% Tested : Matlab 5.1
%     By : Eran O. Ofek                    Nov 1996
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=(1:1:100).'; Y=sin(2.*pi.*X./5)+randn(size(X));
%          S=timeseries.period_entropy([X,Y],[0 1 0.0025],10);
% Reliable: 2
%--------------------------------------------------------------------------

c_x = 1;
c_y = 2;

NumObs = length(x(:,c_x));

LowFreq = FreqVec(1);
HighFreq = FreqVec(2);
FreqInt = FreqVec(3);

% find maximun and minimum of magbnitude for normalization
MinMag = min(x(:,c_y));
MaxMag = max(x(:,c_y));
DelMag = MaxMag - MinMag;

EntMat = zeros(m,m);
S = zeros(length([LowFreq:FreqInt:HighFreq]),2);

N = 1;
for Freq=LowFreq:FreqInt:HighFreq,
   temp1 = x(:,c_x).*Freq;
   phase = temp1 - floor(temp1);
   % calculating the density matrix
   x_pos = floor(phase.*m) + 1;
   y_pos = floor(m.*(x(:,c_y) - MinMag)./(DelMag+1e-5)) + 1;
   % calculating the entropy for the trail period
   for k=1:1:NumObs,
      EntMat(y_pos(k), x_pos(k)) = EntMat(y_pos(k), x_pos(k)) + 1;
   end
   % shold be equal to NumObs
   if sum(sum(EntMat))~=NumObs,
    error('ne');
   end
   mu = EntMat./NumObs;
   S(N,2) = -sum(sum(mu.*log(mu+eps)));
   S(N,1) = Freq;
   N = N + 1;
   EntMat = zeros(m,m);
end

% normalize the Entropy
S(:,2) = S(:,2)./log(m.*m);


