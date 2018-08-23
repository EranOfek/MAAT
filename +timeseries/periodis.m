function [Pxw, Fm] = periodis(x,l_f,h_f,df,pr)
% Scargle periodogram. OBSOLETE: Use period.m instead.
% Package: timeseries
% Description: Lomb-Scargle periodigram. Calculate the Lomb-Scargle power
%              spectrum for a time series.
%              OBSOLETE: Use timeseries.period instead.
% Input  : - Timeseries matrix, [Time, Mag].
%          - Lowest frequency (h_l).
%          - Highest frequency (h_f).
%          - Frequency interval (df).
%          - The probability cutoff (pr), default no cutoff.
%            Power spectra peaks with probability smaller than this
%            cutoff are not listed.
% Output : - Periodigram matrix [Frequency, Power] normalized by the
%            variance of the data points (Horne & Baliunas 1986).
%          - The peaks of the periodogram sorted by the spectral power.
%            [Frequency, Power, Period, (1-False alarm probability)].
% See Also : periodia; periodit; pdm; fitharmo; period
% Reference: Scargle, J.D. ApJ 263, 835-853 (1982).
%	     Horne, J.H. & Baliunas, S.L. ApJ 302, 757-763 (1986).
% Tested : Matlab 4.2
%     By : Eran O. Ofek                    Mar 1994
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 1
%--------------------------------------------------------------------------

c_x = 1;
c_y = 2;
if (nargin==1)
   % set by default
   Range = range(x(:,c_x));
   l_f = 2./Range;
   h_f = 0.5./min(diff(sort(x(:,c_x))));
   df  = 0.25./Range;
elseif nargin==4
   pr = 0; 
elseif nargin==5
   % do nothing
else
   error('Illegal number of input arguments');
end


noise = std(x(:,c_y)).*std(x(:,c_y));
N0 = length(x(:,c_x));
Ni = -6.362 + 1.193.*N0 + 0.00098.*N0.*N0;
tem0 = x(:,c_y) - mean(x(:,c_y));
f_ind = l_f;
F   = [l_f:df:h_f].';
Nf  = length(F);
%Pxw = zeros((h_f-l_f)./df,2);
Pxw = zeros(Nf,2);
Pxw(:,1) = F;

for k=1:1:Nf
   om = 2.*pi.*F(k,1);
   tau = sum(sin(2.*om.*x(:,c_x)))./sum(cos(2.*om.*x(:,c_x)));
   tau = atan(tau)./(2.*om);
   Axc = cos(om.*(x(:,c_x) - tau));
   Axs = sin(om.*(x(:,c_x) - tau));
   Ax1 = sum(tem0.*Axc);
   Ax2 = sum(tem0.*Axs);
   Pxw(k,2) = 0.5.*(Ax1.*Ax1./sum(Axc.*Axc) + Ax2.*Ax2./sum(Axs.*Axs));

end
% normalization of the periodogram
Pxw(:,2) = Pxw(:,2)./noise;



if (nargout==2)
   DiffVec = diff(sign(diff([0;Pxw(:,2);0])));
   IdV     = find(DiffVec==-2);

   Fm      = [Pxw(IdV,1), Pxw(IdV,2), 1./Pxw(IdV,1), (1-exp(-Pxw(IdV,2))).^Ni];
   Fm      = sortrows(Fm,2);
end
