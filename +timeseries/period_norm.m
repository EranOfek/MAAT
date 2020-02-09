function PS=period_norm(Data,FreqVec,Norm,SubMean)
% Normalized power spectrum of non equally spaced time series
% Package: timeseries
% Description: Calculate the normalized normal power spectrum of a times
%              series. See period.m for a more flexiable function.
% Input  : - Two column matrix containing the time series
%            [Time, measurment] or [Time, measurment, error].
%          - This is a column vector of frequencies at which to
%            calculate the power spectrum.
%          - Normalization method:
%            'Var' - Normalize by variance (Default).
%            'Amp' - Normalize by amplitude.
%          - Subtract mean {true|false}. Default is true.
% Output : - Two columns matrix of the un-normalized power spectrum
%            [frequency, power_spec].
% See also: period.m
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    May 2011
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
   PS(FreqInd,2) = abs(sum(M.*exp(-2.*pi.*1i.*T.*FreqVec(FreqInd)))).^2./N;
end

switch lower(Norm)
 case 'amp'
    % do nothing
 case 'no'
    PS(:,2) = PS(:,2).*N;
 case 'var'
    PS(:,2) = PS(:,2)./var(M);
 otherwise
    error('Unknwon normalization option');
end
