function [ShiftedImage,NY,NX,Nr,Nc]=imagefft_shift_fft(ImageFFT,DX,DY,NY,NX,Nr,Nc)
% Shift the fft(Image) using the sub pixel Fourier shift theorem
% Package: ImUtil.Im
% Description: Shift an fft2(image) using the FFT shift thorem. This works
%              well when the image does not contain sharp artifacts.
%              Sharp artifacts will produce ringing.
%              Note that the shift is defined on the content of the image,
%              rather than the image boundries - e.g., the stars will be
%              shifted in the requested direction.
% Input  : - An fft2 of an image (2D matrix).
%          - X shift to apply to input image.
%          - Y shift to apply to input image.
%          - NY (supply for faster performences). See output.
%          - NX (supply for faster performences). See output.
%          - Nr (supply for faster performences). See output.
%          - Nc (supply for faster performences). See output.
% Output : - Shifted image (in real space - not ffted) with the same
%            size as the input image.
%          - NY
%          - NX
%          - Nr
%          - Nc
% See also: SIM/image_shift_fft.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ShiftedImage=ImUtil.Im.image_shift_fft(Image,1.22,-3.1);
% Reliable: 2
%--------------------------------------------------------------------------

Phase = 2;

Nim = 1;
if (nargin<4),
    % NY, NX, Nr, Nc are not provided by user
    [NY,NX,Nim] = size(ImageFFT);
    Nr = ifftshift((-fix(NY/2):ceil(NY/2)-1));
    Nc = ifftshift((-fix(NX/2):ceil(NX/2)-1));
    [Nc,Nr] = meshgrid(Nc,Nr);
end

% Fourier Transform shift theorem
if (Nim==1),
    ShiftedImage = ifft2(ImageFFT.*exp(-1i.*2.*pi.*(DY.*Nr./NY + DX.*Nc./NX))).*exp(-1i.*Phase);
else
    % Image is cube
    % not operational
    %ShiftedImage = ifft2(fft2(Image).*exp(-1i.*2.*pi.*(  bsxfun(@times,DY,shiftdim(Nr,-1))./NY + bsxfun(@times,DX,shiftdim(Nc,-1))./NX))).*exp(-1i.*Phase);
end
ShiftedImage = abs(ShiftedImage);









