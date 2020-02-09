function PS=period_normnl(Data,FreqVec,Norm,SubMean)
% Normzlied power spectrum using no loops (may be faster in some cases)
% Package: timeseries
% Description: Calculate the classical (Lomb) power spectrum of a time
%              series using no loops.
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
% Output : - Two columns matrix of the un-normalized power spectrum
%            [frequency, power_spec].
% See also: period.m
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    May 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: FreqVec=(0:0.001:1).'; T=(1:1:500).'; Data = [T,sin(2.*pi.*0.1.*T)];
%          PS = timeseries.period_normnl(Data,FreqVec);
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

if (SubMean)
    M       = Data(:,Col.M) - mean(Data(:,Col.M));
else
    M       = Data(:,Col.M);
end

% 
%PS = [FreqVec, abs(sum(bsxfun(@times,M.',exp(-2*pi*1i.*FreqVec * T.')),2)).^2./N];         
PS = [FreqVec, abs(sum(bsxfun(@times,M.',exp(-2*pi*1i.*bsxfun(@times,FreqVec,T.'))),2)).^2./N];   % faster!

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
