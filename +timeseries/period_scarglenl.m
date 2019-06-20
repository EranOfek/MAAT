function PS=period_scarglenl(Data,FreqVec,Norm,SubMean)
% Scargle power spectrum of non equally spaced time series / no loops
% Package: timeseries
% Description: Calculate the un-normalized Scargle power spectrum of a times
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
% Reference: Scargle, J.D. ApJ 263, 835-853 (1982).
%            Horne, J.H. & Baliunas, S.L. ApJ 302, 757-763 (1986).
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

M       = Data(:,Col.M) - mean(Data(:,Col.M));
PS      = zeros(Nf,2);
PS(:,1) = FreqVec;
for FreqInd=1:1:Nf
   Tau           = atan(sum(sin(4.*pi.*FreqVec(FreqInd).*T))./sum(cos(4.*pi.*FreqVec(FreqInd).*T)))./(4.*pi.*FreqVec(FreqInd));

   PS(FreqInd,2) = abs(sum(M.*exp(-2.*pi.*1i.*(T-Tau).*FreqVec(FreqInd)))).^2./N;
end

switch lower(Norm)
 case 'amp'
    % do nothing
 case 'var'
    PS(:,2) = PS(:,2)./var(M);
 otherwise
    error('Unknwon normalization option');
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
%PS = [FreqVec, abs(sum(bsxfun(@times,M.',exp(-2*pi*1i.*bsxfun(@times,FreqVec,T.'))),2)).^2./N];   % faster! not in new matlab versions

error('bug - not equal to period_scargle')
Tau = atan(sum(sin(4.*pi.*FreqVec.*T.'))./sum(cos(4.*pi.*FreqVec.*T.')))./(4.*pi.*FreqVec);

Tmp = [2.*pi.*FreqVec.*(T.'-Tau)].';
PS  = [FreqVec, 0.5.*( sum(M.*cos(Tmp) ).^2./sum( cos(Tmp  ).^2 )  + ...
                       sum(M.*sin(Tmp) ).^2./sum( sin(Tmp  ).^2 )  ).'];
       


switch lower(Norm)
 case 'amp'
    % do nothing
 case 'var'
    PS(:,2) = PS(:,2)./var(M);
 otherwise
    error('Unknwon normalization option');
end
