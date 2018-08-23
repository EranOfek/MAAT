function Lap=imlaplacian(Matrix)
% Laplacian filter for a 2D matrix
% Package: ImUtil.Im
% Description: Calculate the laplacian of a 2-D image using a convolution
%              kernel.
%              This function is considerably faster than del2.m
% Input  : - 2-D matrix.
% Output : - The laplacian of the image.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: R=randn(2048,2048)+100; R(900,700)=1000;
%          Lap=ImUtil.Im.imlaplacian(R)
% Reliable: 2
%--------------------------------------------------------------------------

LaplaceKernel = [0 -1 0; -1 4 -1; 0 -1 0];
Size = size(Matrix);
% will be faster to use filter2...
%Lap  = (abs(ifft2(fft2(Matrix).*fft2(LaplaceKernel,Size(1),Size(2)))));
Lap = filter2(LaplaceKernel,Matrix);