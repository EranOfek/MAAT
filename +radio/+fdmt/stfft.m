function S=stfft(raw_signal,block_size)
% Short time fast fourier transform
% Package: radio.fdmt
% Description: Short time fast fourier transform (STFFT) of a time series.
%              Each block (time bin) in the input series is FFTed,
%              and the output is the FFT in each block. For a raw radio
%              signal of voltage as a function of time, this function
%              will return the amplitude as a function of frequency in each
%              time bin.
% Input  : - The raw signal series. A vector of measurments.
%          - Block size - number of measurments in each block.
% Output : - A matrix of time vs. frequency.
%            Note that absolute value squared is not performed.
% Tested : Matlab R2014a
%     By : Barak Zackay                    Nov 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=radio.fdmt.stfft(rand(10000,1),100);
% Reliable: 2
%--------------------------------------------------------------------------

S = raw_signal(1:floor(length(raw_signal)/block_size) * block_size );
S = reshape(S,[block_size floor(length(raw_signal)/block_size)]);
S = fft(S,[],1);
S = transpose(S);
