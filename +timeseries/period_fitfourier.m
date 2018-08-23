function [PS,Par,Stat]=period_fitfourier(Data,FreqVec,Harmon,PolyN)
% Fit a Fourier series to a time series
% Package: timeseries
% Description: Fit a polynomial and fourier series as a function of
%              frequency to a time series.
%              The fitted series is:
%              Y = A_ + A_*T + A_*T^3,... A_*sin(2*pi*f*T*H(1)) +
%                                         A_*sin(2*pi*f*T*H(2)) + ...
%                                         A_*cos(2*pi*f*T*H(1)) +
%                                         A_*cos(2*pi*f*T*H(2)) + ...
%              See period.m for a more flexiable function.
% Input  : - Two column matrix containing the time series
%            [Time, measurment] or [Time, measurment, error].
%          - Frequency range and interval in which to calculate the
%            power spectrum.
%            This is a column vector of frequencies at which to
%            calculate the power spectrum.
%          - Row vector of Harmonies to fit, e.g. [1 2 3].
%            Default is [1 2].
%          - Degree of polynomials to fit. Default is 1.
% Output : - Two columns matrix of the power spectrum equivalent
%            [frequency, amplitude], where amplitude is the total amplitude
%            of all the harmonies at a given frequency.
%          - Structure containing the best fit parameters and errors
%            for each frequency.
%            [1 T T^2 T^3,... sin(2*pi*f*H*T), ..., cos(2*pi*f*H*T),...]
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    May 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [PS,Par,Stat]=timeseries.period_fitfourier(Data,FreqVec,Harmon,PolyN)
%----------------------------------------------------------------------------

Def.Harmon  = [1 2];
Def.PolyN   = 1;
if (nargin==2),
   Harmon = Def.Harmon;
   PolyN  = Def.PolyN;
elseif (nargin==3),
   PolyN  = Def.PolyN;
elseif (nargin==4),
   % do nothing
else
    error('Illegal number of input arguments');
end

Col.T = 1;
Col.M = 2;
Col.E = 3;

T       = Data(:,Col.T);
N       = numel(T);
Nf      = numel(FreqVec);

% construct the design matrix:
Hpoly = [];
for Ip=0:1:PolyN,
   Hpoly = [Hpoly, T.^Ip];
end
if (size(Data,2)==2),
   Data = [Data, ones(N,1)];
end

InvVar = 1./Data(:,Col.E).^2;   % inverse variance

Stat.Npar = PolyN+1+length(Harmon).*2;
Stat.Dof  = N - Stat.Npar;

% fit harmonies for each frequency
Chi2      = zeros(Nf,1);
RMS       = zeros(Nf,1);
Par.Par    = zeros(Nf,Stat.Npar);
Par.Err = zeros(Nf,Stat.Npar);

for FreqInd=1:1:Nf,
   Hharm = [sin(2.*pi.*FreqVec(FreqInd).*(T*Harmon)),...
            cos(2.*pi.*FreqVec(FreqInd).*(T*Harmon))];

   H = [Hpoly, Hharm];
   [P,E] = lscov(H,Data(:,Col.M),InvVar);
   Resid = Data(:,2) - H*P;
   Stat.Chi2(FreqInd) = sum((Resid./Data(:,3)).^2);
   Stat.RMS(FreqInd) = std(Resid);

   Par.Par(FreqInd,:)    = P.';
   Par.Err(FreqInd,:) = E.';
end

% amplitude normalization
PS = sqrt(sum(Par.Par(:,PolyN+1:end).^2,2));

