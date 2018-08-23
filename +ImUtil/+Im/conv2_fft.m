function Res=conv2_fft(Mat1,Mat2)
%--------------------------------------------------------------------------
% conv2_fft function                                               ImBasic
% Description: Convolve two equal size matrices using fft.
%              This is somewhat faster than conv_fft2.m.
% Input  : - A matrix.
%          - A matrix the same size  as the first input.
% Output : - Convolution between the two matrices.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: conv2_fft(rand(1000,1000),rand(1000,1000));
% Reliable: 
%--------------------------------------------------------------------------

Size1 = size(Mat1);
Size2 = size(Mat2);

if (any(Size2>Size1))
    error('Second matrix must be equal or smaller than first matrix - use conv_fft2.m instead');
end

% padding
MaxSize = max(Size1,Size2);
MinSize = min(Size1,Size2);
Size    = MaxSize + floor(MinSize.*0.5);
Mat1    = padarray(Mat1,Size-Size1,'post');
Mat2    = padarray(Mat2,Size-Size2,'post');
Sh1     = floor(MinSize.*0.5);

% convoltion
Res = ifft2(fft2(Mat1).*fft2(Mat2));      

% remove padding
Res = Res(1+Sh1(1):end,1+Sh1(2):end);

% if (Size1(1)<Size2(1)),
%      Res = Res(1:end-1,:);
% elseif (Size1(1)>Size2(1)),
%      Res = Res(2:end,:);
% end
% if (Size1(2)<Size2(2)),
%      Res = Res(:,1:end-1);
% elseif (Size1(2)>Size2(2)),
%      Res = Res(:,2:end);    
% end


