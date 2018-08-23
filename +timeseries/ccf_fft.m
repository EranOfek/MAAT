function CCF=ccf_fft(X1,X2)
% Cross correlation function of evenly spaced data using fft
% Package: timeseries
% Description: Calculate the normalized cross correlation function of
%              two evenly spaced timeseries using fft.
% Input  : - First series. Either Y or [X,Y].
%          - Second series. Either Y or [X,Y].
% Output : - Structure containing the following fields:
%            .X - The lag
%            .C - The normalized cross correlation.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Dec 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: R=randn(100,1); CCF=timeseries.ccf_fft(R,R); %plot(CCF.X,CCF.C)
%          %An example for a matched filter use:
%          Std = 0.5;
%          R=randn(100000,1).*Std; % generate a random equally spaced time series.
%          T=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1]'; % define a filter
%          R(30001:30015)=R(30001:30015) + T; % insert the template to the data
%          tic;CCF=timeseries.ccf_fft(R,T);toc
%          %plot(CCF.X,CCF.C); % you can clearly see the peak at ~30,000
% Reliable: 2
%--------------------------------------------------------------------------
Col.X = 1;
Col.Y = 2;
[N1,C1] = size(X1);
[N2,C2] = size(X2);

if (C1==1)
   X1 = [(1:1:N1).', X1];
end
if (C2==1)
   X2 = [(1:1:N2).', X2];
end
D = diff(X1(:,Col.X));
D = D(1);

N = max(N1,N2);

   

Norm = N.*std(X1(:,Col.Y)).*std(X2(:,Col.Y));
CCF.C = (ifft(fft(X1(:,Col.Y)-mean(X1(:,Col.Y)),N).*conj(fft(X2(:,Col.Y)-mean(X2(:,Col.Y)),N)))) ./Norm;
CCF.X = (0:1:N-1).'.*D;

