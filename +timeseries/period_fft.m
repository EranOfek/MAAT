function [PS,Peaks]=period_fft(Data,Normalization)
% Power spectrum of evenly spaced time series using fft
% Package: timeseries
% Description: Calculate power spectrum for evenly spaced time series
%              using fast Fourier transform.
%              See also period.m
% Input  : - Two columns matrix (e.g., [time, x]).
%          - Power spectrum normalization:
%            'no'    - no normalization.
%            'amp'   - amplitude normalization (e.g., sqrt of power).
%            'var'   - Normalize by variance, default.
% Output : - Power spectrum [Frequency, power].
%          - Structure containing information about the peaks in the power
%            spectrum. The following fields are available:
%            .FAB   - 1 - False alarma probability
%            .Freq  - Frequency
%            .Per   - Period (1/Frequency)
%            .PS    - Power.
% See also: period.m
% Tested : Matlab 7.6
%     By : Eran O. Ofek                    Jul 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------
Def.Normalization = 'var';
if (nargin==1),
   Normalization = Def.Normalization;
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end

DT = Data(2,1) - Data(1,1);  % Time diff. between measurments

Col.T = 1;
Col.Y = 2;
N = size(Data,1);
FFT = fft(Data(:,Col.Y));

Freq = ((1:1:N).'-1)./N;
Freq = Freq./DT;

PS  = [Freq, abs(FFT).^2];

switch lower(Normalization)
 case 'no'
    % no normalization
    % do nothing
 case 'amp'
    % amplitude normalization
    PS = [PS(:,1), sqrt(PS(:,2))./N];
 case 'var'
    % variance normalization
    Var = std(Data(:,Col.Y)).^2;
    PS = [PS(:,1), PS(:,2)./Var./N];
 otherwise
    error('Unknown Normalization option');
end



if (nargout>1),
   Ipeak = diff(sign(diff([0;PS(:,2);0])))==-2;
   Nd    = size(Data,1);
   Ni = -6.362 + 1.193.*Nd + 0.00098.*Nd.^2;  % number of independent data points : Horne, J.H. & Baliunas, S.L. ApJ 302, 757-763 (1986)

   Peaks.FAB  = (1-exp(-PS(Ipeak,2))).^Ni;
   %[Peaks.FAB, SI] = sort(Peaks.FAB);
   Peaks.Freq = PS(Ipeak,1);
   Peaks.Per  = 1./PS(Ipeak,1);
   Peaks.PS   = PS(Ipeak,2);
   Sorted = sortrows([Peaks.FAB, Peaks.Freq, Peaks.Per, Peaks.PS],1);
   Peaks.FAB  = Sorted(:,1);
   Peaks.Freq = Sorted(:,2);
   Peaks.Per  = Sorted(:,3);
   Peaks.PS   = Sorted(:,4);
  
end

