function [PS,FreqVec]=period_norm_bootstrap(Data,FreqVec,Norm,SubMean,Nsim)
% Run bootstrap normal pertiodogram on a light curve.
% Package: timeseries
% Description: 
% Input  : - Two column matrix containing the time series
%            [Time, measurment] or [Time, measurment, error].
%          - This is a column vector of frequencies at which to
%            calculate the power spectrum.
%          - Normalization method:
%            'Var' - Normalize by variance (Default).
%            'Amp' - Normalize by amplitude.
%          - Subtract mean {true|false}. Default is true.
%          - Number of bootstrap simulations. Default is 1000.
% Output : - Matrix in wgich each column is the power spectrum for a
%            specific bootstrap simulation. Number of columns is the number
%            of simulations.
%          - Vector of corresponding frequencies.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: t=timeseries.random_time_sequence;
%          x=0.2.*sin(2.*pi.*t./11) + randn(numel(t),1);
%          F = (0:0.001:1)';
%          P=timeseries.period([t,x],F);
%          [PS,F]=timeseries.period_norm_bootstrap([t,x],F);
%          plot(F,P(:,2)./std(PS,[],2))
% Reliable: 2
%--------------------------------------------------------------------------



%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Def.Norm    = 'Var';
Def.SubMean = true;
Def.Nsim    = 1000;
if (nargin==2)
    Norm    = Def.Norm;
    SubMean = Def.SubMean;
    Nsim    = Def.Nsim;
elseif (nargin==3)
    SubMean = Def.SubMean;
    Nsim    = Def.Nsim;
elseif (nargin==4)
    Nsim    = Def.Nsim;
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

% bootstrap M
Iboot = randi(N,[N, Nsim]);
Mboot = M(Iboot);

PS      = zeros(Nf,Nsim);
PS(:,1) = FreqVec;
for FreqInd=1:1:Nf
   PS(FreqInd,:) = abs(sum(Mboot.*exp(-2.*pi.*1i.*T.*FreqVec(FreqInd)),1)).^2./N;
end

switch lower(Norm)
 case 'amp'
    % do nothing
 case 'var'
    PS(:,2) = PS(:,2)./var(M);
 otherwise
    error('Unknwon normalization option');
end
