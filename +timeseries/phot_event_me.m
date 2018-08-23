function [S,Sig]=phot_event_me(x,FreqVec,m)
%--------------------------------------------------------------------
% phot_event_me function     Searching periodicity in time-tagged
%                     events using information entropy.
%                     For each trail period, the phase-magnitude space
%                     is divided into m by m cells, and the information
%                     entropy is calculated.
%                     Then the probability to find observations
%                     in each square is MU(i) and the Entropy is
%                     (-Sigma{MU*ln(MU)}).
%                     The output Entropy is normalized by ln(m*m).
% input  : - Sorted vector of events.
%          - Frequency Vector, [Lower Freq., Higher Freq., Freq. Interval]
%          - number of elements in phase (number of bins). 
% output : - Entropy vs. Frequency.
%            This is a two column matrix. The first column contains
%            the frequency and the second list the entropy.
%          - One sigma probability.
% Reference : Cincotta et al. 1999, MNRAS 302,582.
% Tested : Matlab 5.2
%    By  Eran O. Ofek           April 1999
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------


NumObs = length(x);

LowFreq = FreqVec(1);
HighFreq = FreqVec(2);
FreqInt = FreqVec(3);


EntVec = zeros(m);
S = zeros(length([LowFreq:FreqInt:HighFreq]),2);

N = 1;
for Freq=LowFreq:FreqInt:HighFreq
   phase = x.*Freq - floor(x.*Freq);
   % calculating the density vector
   x_pos = floor(phase.*m) + 1;
   % calculating the entropy for the trail period
   for k=1:1:NumObs
      EntVec(x_pos(k)) = EntVec(x_pos(k)) + 1;
   end
   % shold be equal to NumObs
   if sum(sum(EntVec))~=NumObs
    error('ne');
   end
   mu = EntVec./NumObs;
   S(N,2) = -sum(sum(mu.*log(mu+eps)));
   S(N,1) = Freq;
   N = N + 1;
   EntVec = zeros(m);
end

% normalize the Entropy
S(:,2) = S(:,2)./log(m);

% statistics - estimate one sigma
Sig = sqrt(1./NumObs).*((log(NumObs)+1)./log(m) - 1);


