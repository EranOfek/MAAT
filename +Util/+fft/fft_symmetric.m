function [newFFT,Freq]=fft_symmetric(fftD,Dim)
% Make a 1-D fft a complex-conjugate symmetric
% Package: Util.fit
% Description: Given a 1-D fft of a vector, or a matrix, force the negative
%              frequency signal to be complex conjugate of the positive
%              frequency signal.
% Input  : - fft of a data - e.g., fft(data).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [newFFT,Freq]=Util.fft.fft_symmetric(fftD,Dim)
% Reliable: 2
%--------------------------------------------------------------------------



if nargin<2
    Dim = 1;
end
[Ni,Nj] = size(fftD);
if (Ni==1 || Nj==1)
    fftD = fftD(:);
end

N       = size(fftD,Dim);



Freq    = Util.fft.fft_freq(N)
Nf      = numel(Freq);

IF0     = find(Freq==0);


NewInd = [(1:IF0-1), IF0, (IF0-1:-1:1)];
NewInd = NewInd(1:Nf);




if (Dim==1)
    newFFT  = fftD(NewInd,:);
elseif (Dim==2)
    newFFT  = fftD(:,NewInd);
else
    error('Suported Dim are 1 or 2 only');
end


