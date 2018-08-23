function PS=period_norm_order2(Data,FreqVec,Norm,SubMean)
% 2nd order harmonic power spectrum of non evenly spaced time series
% Package: timeseries
% Description: Calculate the normalized power spectrum of the second order
%              of a times series (i.e., assuming that the first and second
%              harmonies have the same amplitude).
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
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    May 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X = (0:1:100).'; Y = sin(2.*pi.*X./11) + sin(4.*pi.*X./11 + 1) +randn(size(X)).*0.01;
%          PS=timeseries.period_norm_order2([X,Y],(0:0.001:1).');
% Reliable: 2
%--------------------------------------------------------------------------
W1 = 1;
W2 = 1;

Def.Norm    = 'Var';
Def.SubMean = true;
if (nargin==2),
    Norm    = Def.Norm;
    SubMean = Def.SubMean;
elseif (nargin==3),
    SubMean = Def.SubMean;
else
    % do nothing
end

Col.T = 1;
Col.M = 2;
T       = Data(:,Col.T);
N       = numel(T);
Nf      = numel(FreqVec);

if (SubMean),
    M       = Data(:,Col.M) - mean(Data(:,Col.M));
else
    M       = Data(:,Col.M);
end
PS      = zeros(Nf,2);
PS(:,1) = FreqVec;
for FreqInd=1:1:Nf,
   PS(FreqInd,2) = abs(W1.*sum(M.*exp(-2.*pi.*1i.*T.*FreqVec(FreqInd))) + ...
                       W2.*sum(M.*exp(-4.*pi.*1i.*T.*FreqVec(FreqInd))) ).^2./N;
end

switch lower(Norm)
 case 'amp'
    % do nothing
 case 'var'
    PS(:,2) = PS(:,2)./var(M);
 otherwise
    error('Unknwon normalization option');
end
