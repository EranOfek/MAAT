function Out=conv2_cell(Image,Kernel,Method)
% Convolve images in a cell array with kernels in a cell array.
% Package: @ImUtil.Im
% Description: Convolve each image in a cell array of images with a
%              corresponding kernel, in a cell array of kernels.
% Input  : - A cell array of images, or a single image (2D array).
%          - A cell array of kernels, or a single kernel (2D array).
%          - Convolution method:
%            'auto'     - Guess the fastest option. Default. 
%            'conv2'    - Use conv2.m
%            'conv2_fft'- Use ImUtil.Im.conv2_fft.m
%            'conv_fft2'- Use Util.external.conv_fft2.m
% Output : - A cell array of convolved images.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Out=conv2_cell(rand(10,10),rand(10,10));
%          Out=conv2_cell({rand(10,10),rand(11,11)},rand(10,10));
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==2)
    Method  = 'auto';
end

if (~iscell(Image))
    Image = {Image};
end
if (~iscell(Kernel))
    Kernel = {Kernel};
end

switch lower(Method)
    case 'conv2'
        FunConv = @conv2;
        Pars    = {'same'};
    case 'conv_fft2'
        FunConv = @Util.external.conv_fft2;
        Pars    = {'same'};
    case 'conv2_fft'
        FunConv = @ImUtil.Im.conv2_fft;
        Pars    = {};
    case 'auto'
        if (max(size(Image{1})./size(Kernel{1}))>10)
            % use conv2.m
            FunConv = @conv2;
            Pars    = {'same'};
        else
            if (numel(Image{1})<1000 && numel(Kernel{1})<1000)
                % use conv2.m as images are small stamps
                FunConv = @conv2;
                Pars    = {'same'};
            else
                FunConv = @ImUtil.Im.conv2_fft;
                Pars    = {};
            end
        end
    otherwise
        error('Unknown Method option');
end

     
Nim = numel(Image);
Nk  = numel(Kernel);
N   = max(Nim,Nk);

Out = cell(size(Image));
for I=1:1:N
    Iim = min(I,Nim);
    Ik  = min(I,Nk);
    
    Out{I} = FunConv(Image{Iim},Kernel{Ik},Pars{:});
end
