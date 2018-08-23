function PS=period_complex(Data,FreqVec,Norm,SubMean)
% Fourier transform of non equally spaced time series
% Package: timeseries
% Description: Calculate the complex fourier transform of an unevenly
%              spaced time series.
%              See period.m for a more flexiable function.
% Input  : - Two column matrix containing the time series
%            [Time, measurment] or [Time, measurment, error].
%          - Frequency range and interval in which to calculate the
%            power spectrum.
%            This is a column vector of frequencies at which to
%            calculate the power spectrum.
%          - Normalization method:
%            'Var' - Normalize by variance (Default).
%            'Amp' - Normalize by amplitude.
%          - Subtract mean {true|false}. Default is true.
% Output : - Two columns matrix of the un-normalized power spectrum
%            [frequency, power_spec].
% See also: period.m
% Tested : Matlab 2015b
%     By : Eran O. Ofek                    Dec 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

Def.Norm    = 'Var';
Def.SubMean = true;
if (nargin==2)
    Norm    = Def.Norm;
    SubMean = Def.SubMean;
elseif (nargin==3)
    SubMean = Def.SubMean;
else
    % do nothing
end

Col.T = 1;
Col.M = 2;
T       = Data(:,Col.T);
N       = numel(T);
Nf      = numel(FreqVec);

if (SubMean)
    M       = Data(:,Col.M) - mean(Data(:,Col.M));
else
    M       = Data(:,Col.M);
end
PS      = zeros(Nf,2);
PS(:,1) = FreqVec;
for FreqInd=1:1:Nf
   PS(FreqInd,2) = sum(M.*exp(-2.*pi.*1i.*T.*FreqVec(FreqInd)))./N;
end

switch lower(Norm)
 case 'amp'
    % do nothing
 case 'var'
    PS(:,2) = PS(:,2)./var(M);
 otherwise
    error('Unknwon normalization option');
end
