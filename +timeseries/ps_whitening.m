function [WTS]=ps_whitening(TS)
%--------------------------------------------------------------------------
% ps_whitening function                                         timeseries
% Description: Whiten a time series using its power spectrum.
%              I.e., ifft(fft(x)./abs(fft(x)))/
% Input  : - Time series [T, X], or [X] if the time series is equally
%            spaced.
% Output : - Whitned time series.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [WTS]=ps_whitening(TS
% Reliable: 
%--------------------------------------------------------------------------




Ncol  = size(TS,2); % number of columns

if (Ncol==1),
    % assume the time series is equally spaced
    FT  = fft(TS);               % fft of time series
    PS  = FT.*conj(FT);          % abs val square - power spectrum
    WTS = ifft(FT./sqrt(PS));    % whitened and ifft
    
else
    % assume the time series is NOT equally spaced
    error('not equally spaced time series are not supported yet');
    
end



