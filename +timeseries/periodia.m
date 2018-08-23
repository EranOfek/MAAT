function [Pxw,Fm]=periodia(x,l_f,h_f,df,c_x,c_y)
% Classical power spectrum of non-evenly space time series (OBSOLETE)
% Package: timeseries
% Description: Classical power-spectrum of a non-evenly spaced time series.
%              The power spectrum is normalized by the variance of the data.
%              OBSOLETE: Use timeseries.period instead.
% Input  : - Time series matrix, [Time, Mag], in which the first column
%            is the time and the second column is the magnitude.
%          - Lowest frequency (h_l).
%          - Highest frequency (h_f).
%          - Frequency interval (df).
%          - The column of time (c_x), default is 1.
%          - The column of magnitudes (c_y), default is 2.
% Output : - Periodigram matrix. normalized with the variance 
%            of the observations (Horne & Baliunas 1986).
%            The first column is the frequency and
%            the second column is the power.
%          - The peaks of the periodogram sorted by the probability.
%            [Frequency, Power, Period, (1-False alarm probability)].
%            The probability is good only for evenly spaced data
% See Also : periodis; periodit; pdm; fitharmo; period
% Reference: Koen, C. 1990, ApJ 348, 700-702.
% Tested : Matlab 4.2
%     By : Eran O. Ofek                    Dec 1993
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 1
%--------------------------------------------------------------------------

if nargin==4,
   c_x = 1;
   c_y = 2;
elseif nargin==5,
   c_y = 2;
elseif nargin==6,
   % do nothing
else
   error('Illegal number of input arguments');
end

% set probability cutoff to zero.
pr = 0;
N0 = length(x(:,c_x));
Ni = -6.362 + 1.193.*N0 + 0.00098.*N0.*N0;
tem0 = x(:,c_y) - mean(x(:,c_y));
f_ind = l_f;
k = 1;
Pxw = zeros((h_f-l_f)./df,2);
while f_ind<h_f,
   Pxw(k,1) = f_ind; 
   temp = abs(sum(tem0.*exp(-i*2*pi*f_ind*x(:,c_x))));
   Pxw(k,2) = temp*temp/N0;
   f_ind = f_ind + df;
   k = k + 1;
end

% normalization of the periodogram
noise = std(x(:,c_y)).^2;
Pxw(:,2) = Pxw(:,2)./noise;


if (nargout==2),
   DiffVec = diff(sign(diff([0;Pxw(:,2);0])));
   IdV     = find(DiffVec==-2);

   Fm      = [Pxw(IdV,1), Pxw(IdV,2), 1./Pxw(IdV,1), (1-exp(-Pxw(IdV,2))).^Ni];
   Fm      = sortrows(Fm,2);
end
