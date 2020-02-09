function FreqVec=fft_freq(N)
% Return the frequencies corresponding to fftshift(fft(vec_of_size_N))
% Package: Util.fft
% Description: Return the frequencies corresponding to
%              fftshift(fft(vec_of_size_N)), without dviding by the total
%              time span.
% Input  : - Number of points.
% Output : - Vector of frequencies.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: FreqVec=Util.fft.fft_freq(4)
% Reliable: 2
%--------------------------------------------------------------------------


if (N.*0.5)==floor(0.5.*N)
    % even
    FreqVec = (-N.*0.5:1:N.*0.5-1).';
else
    FreqVec = (-N.*0.5+0.5:1:N.*0.5-0.5).';
end

