function Res=filter2_fft(Mat1,Mat2)
% Filter an image using fft.
% Package: Util.filter
% Description: Filter an image using fft.
%              Note that this is like the built in filter2.m function but
%              uses fft. This function should be faster than filter2.m
%              when the filter size is larger than ~5% in each dimension
%              from the image.
%              Note that in filter2.m the filter is the first input
%              argument, while here it is the second.
%              See also filter2_fftc.m for even faster version.
% Input  : - A 2-D image (matrix).
%          - A filter (matrix) with the same size as the first input image.
%            The filter is supposed to be centered.
% Output : - A filtered image.
% See also: filter2_fftc.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: filter2_fft(rand(1024,1024),rand(1024,1024));
%          filter2_fft(rand(1024,1024),rand(1024,1024),false);
% Reliable: 2
%--------------------------------------------------------------------------


Mat2 = rot90(Mat2,2);

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

% In order to use conj instead of rot
% Mat1 and Mat2 must be the same size without padding!!
% so for the general case do not use conj...
%Res = ifft2(fft2(Mat1).*conj(fft2(fftshift(Mat2))));      

% remove padding
Res = Res(1+Sh1(1):end,1+Sh1(2):end);
