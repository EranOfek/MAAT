function [ShiftedImage,NY,NX,Nr,Nc]=image_shift_fft(Image,DX,DY,NY,NX,Nr,Nc)
% Shift Image using the sub pixel Fourier shift theorem (sinc interp.)
% Package: ImUtil.Im
% Description: Shift an image using the FFT shift thorem. This works well
%              when the image does not contain sharp artifacts.
%              Sharp artifacts will produce ringing.
%              Note that the shift is defined on the content of the image,
%              rather than the image boundries - e.g., the stars will be
%              shifted in the requested direction.
% Input  : - An image (2D matrix).
%          - X shift to apply to input image.
%          - Y shift to apply to input image.
%          - NY (supply for faster performences). See output.
%          - NX (supply for faster performences). See output.
%          - Nr (supply for faster performences). See output.
%          - Nc (supply for faster performences). See output.
% Output : - Shifted image with the same size as the input image.
%          - NY
%          - NX
%          - Nr
%          - Nc
% See also: ImUtil.Im.imagefft_shift_fft.m, SIM/image_shift_fft.m,
%           SIM/imagefft_shift_fft.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ShiftedImage=ImUtil.Im.image_shift_fft(Image,1.22,-3.1);
% Reliable: 2
%--------------------------------------------------------------------------

NewAlgo = true;

if (NewAlgo)
  
    [NY,NX] = size(Image);
    if (0.5*NX==floor(0.5*NX))
        % NX is even 
        PadX = 0;
    else
        PadX = 1;
        NX   = NX + 1;
    end
    if (0.5*NY==floor(0.5*NY))
        % NX is even 
        PadY = 0;
    else
        PadY = 1;
        NY   = NY + 1;
    end
    Image = padarray(Image,[PadX PadY],0,'post');

    % Kernel for X dimension
    OperX = fft([0 1 zeros(1,NX-2)]);
    KernelX = fftshift(exp(1i.*DX.*phase(OperX)));
    KernelX = KernelX./KernelX(1);
    KernelX(NX.*0.5+1) = 1;
    %KernelX = ifft(KernelX);

    % Kernel for Y dimension
    OperY = fft([0 1 zeros(1,NY-2)]);
    KernelY = fftshift(exp(1i.*DY.*phase(OperY))).';
    KernelY = KernelY./KernelY(1);
    KernelY(NY.*0.5+1) = 1;
    %KernelY = ifft(KernelY);

    
    SX = ifft( bsxfun(@times,fft(Image,[],2),KernelX) ,[],2);
    % need to take the real part as there is some residual imaginary
    % part due to computer precision errors
    ShiftedImage=real(ifft( bsxfun(@times,fft(SX,[],1), KernelY) ,[],1));
    
    %Kernel = KernelY*KernelX;
    %fftKernel = KernelY*KernelX;
    
    %ShiftedImage = real(ifft2(fft2(Image).*fft2(Kernel)));
    %ShiftedImage = real(ifft2(fft2(Image).*fftKernel));
    if (PadX==1)
        ShiftedImage = ShiftedImage(:,1:end-1);
    end
    if (PadY==1)
        ShiftedImage = ShiftedImage(1:end-1,:);
    end
    Nr = [];
    Nc = [];
else
    %function [ShiftedImage,NY,NX,Nr,Nc]=image_shift_fft(Image,DX,DY,NY,NX,Nr,Nc)
    % old algorithm - better
    Phase = 2;

    % add bias to avoid negative numbers
    MinVal = min(Image(:))+1;
    %MinVal = 0;
    Image  = Image + MinVal;

    % [NY1,NX1] = size(Image);
    % Image = padarray(Image,[NY1 NX1],0,'both');
    % [NY,NX,Nim] = size(Image);

    Nim = 1;
    if (nargin<4)
        % NY, NX, Nr, Nc are not provided by user
        [NY,NX,Nim] = size(Image);

        Nr = ifftshift((-fix(NY/2):ceil(NY/2)-1));
        Nc = ifftshift((-fix(NX/2):ceil(NX/2)-1));
        [Nc,Nr] = meshgrid(Nc,Nr);
    end

    % Fourier Transform shift theorem
    if (Nim==1)  
        ShiftedImage = ifft2(fft2(Image).*exp(-1i.*2.*pi.*(DY.*Nr./NY + DX.*Nc./NX))).*exp(-1i.*Phase);
    else
        % Image is cube
        error('Cube images not supported yet');
        %ShiftedImage = ifft2(fft2(Image).*exp(-1i.*2.*pi.*(  bsxfun(@times,DY,shiftdim(Nr,-1))./NY + bsxfun(@times,DX,shiftdim(Nc,-1))./NX))).*exp(-1i.*Phase);
    end

    % add bias value to image
    ShiftedImage = abs(ShiftedImage) - MinVal;
    %ShiftedImage = ShiftedImage(NY1+1:2*NY1, NX1+1:2*NX1);
end
