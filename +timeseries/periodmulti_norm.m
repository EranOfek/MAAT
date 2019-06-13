function [FreqVec,PS]=periodmulti_norm(T,Data,FreqVec,Norm)
% Simultanous power spectrum of time series with common times
% Package: timeseries
% periodmulti_norm function                                     timeseries
% Description: Calculate the normal power spectrum of a set of multiple
%              times series which have common times.
%              This program run the power spectrum simultaneously for all
%              time serieses and its is faster than running the power
%              spectrum on a single time series at a time.
% Input  : - Vector of common times for all time serieses.
%          - Matrix of measurements in which each column refer to one
%            time series and each row refer to one time.
%          - Vector of frequencies in which to calculate the power spectrum.
%          - Normalization method:
%            'Var' - Normalize by variance (Default).
%            'Amp' - Normalize by amplitude.
% Output : - Vector of frequencies.
%          - Matrix of power spectrum. Each column represent the power
%            spectrum for each time series.
% See also: period.m, period_norm.m
% Tested : Matlab 7.16
%     By : Eran O. Ofek                    Sep 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

Def.Norm = 'Var';
if (nargin==3)
   Norm = Def.Norm;
end

%Col.T = 1;
%Col.M = 2;
%T       = Data(:,Col.T);
N       = numel(T);
Nf      = numel(FreqVec);
Ns      = size(Data,2);   % number of columns (i.e., time series)

M       = bsxfun(@minus,Data,mean(Data));
PS      = zeros(Nf,Ns);
%PS(:,1) = FreqVec;
for FreqInd=1:1:Nf
   PS(FreqInd,:) = abs(nansum(bsxfun(@times,M,exp(-2.*pi.*1i.*T.*FreqVec(FreqInd))))).^2./N;
end

switch lower(Norm)
 case 'amp'
    % do nothing
 case 'var'
	   PS = bsxfun(@rdivide,PS,var(M));
 otherwise
    error('Unknwon normalization option');
end
