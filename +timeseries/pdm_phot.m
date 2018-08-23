function [F,Variance,StD_Mean]=pdm_phot(X,Nb,Freq_Vec)
%--------------------------------------------------------------------
% pdm_phot function    Phase Dispersion Minimization to Photon
%                    arrival time data.
% Input  : - Sorted vector of arrival times.
%          - Number of bins.
%          - Vector of frequencies to search,
%            [Low Frequency, High Frequency, Frequency Interval]
%            default is to choose the Low Freq. as 1/(Time Span)
%            the High Freq. as the mean(diff(photon arrival time))
%            Freq. Interval as 1/(4 X Time Span)
% Output : - Vector of frequencies
%          - Vector of variances coresponding to each trial frequency.
%          - Vector of standard deviations of means coresponding
%            to each trial frequency.
% Tested : Matlab 5.1
%     By : Eran O. Ofek           November 1996
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------
if nargin==2,
   Time_Span = Max(X) - Min(X);
   Low_Freq  = 1./(Time_Span);
   High_Freq = mean(diff(X));
   Freq_Int  = 1./(4.*Time_Span);
elseif nargin==3,
   % do nothing
else
   error('Illegal number of input arguments');
end


% assuming X is sorted
%X = sort(X)

% Number of Photons
NoP = length(X);


Low_Freq  = Freq_Vec(1);
High_Freq = Freq_Vec(2);
Freq_Int  = Freq_Vec(3);


% calculating Variance for all the data
Diff_All = diff(X);
Vari_All = sum((Diff_All - mean(Diff_All)).^2)./(NoP - 2);

% initialize Vector of Indices
%VI = zeros(1,Nb-1);
% initialize Vector of StD
StD_V = zeros(1,Nb);
% initialize Vector of Means
Mean_V = zeros(1,Nb);
% initialize Vector of Number of Photons in Bin
NoPB_V = zeros(1,Nb);
% initialize Vector of Variance
Variance = zeros(1,length(Low_Freq:Freq_Int:High_Freq));
% initialize Vector of StD of Means 
StD_Mean = zeros(1,length(Low_Freq:Freq_Int:High_Freq));
% initialize Vector of Frequencies of Means 
F = [Low_Freq:Freq_Int:High_Freq];


J = 1;
for Freq=Low_Freq:Freq_Int:High_Freq,
   F(J)      = Freq;
   P         = X.*Freq - floor(X.*Freq);
   Sorted_P  = sort(P);
   Diff_P    = diff(Sorted_P);

   First_Ind = 1;
   for Bin_N=1:1:Nb-1,
      Bin_Ind       = bin_sear(Sorted_P, Bin_N./Nb);
      StD_V(Bin_N)  = std(First_Ind:Bin_Ind);
      Mean_V(Bin_N) = mean(First_Ind:Bin_Ind);
      NoPB_V(Bin_N) = Bin_Ind -  First_Ind + 1;

      First_Ind     = Bin_Ind
   end
   StD_V(Nb)  = std(First_Ind:NoP);
   Mean_V(Nb) = mean(First_Ind:NoP);
   NoPB_V(Nb) = NoP - First_Ind + 1;


   % calculate statistics in each bin
   Variance(J) = sum((NoPB_V - 1).*StD_V.^2./(sum(NoPB_V) - Nb);
   StD_Mean(J) = std(Mean_V);
end

   



