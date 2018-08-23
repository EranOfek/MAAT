function [Val,Image]=image_posval(Image,X,Y,Filter,FiltFun)
% Get values of image pixels at specific rounded X/Y positions.
% Package: ImUtil.Im
% Description: Read values of image pixels at specific rounded
%              X/Y positions.
%              Optionally, filter the image.
% Input  : - Image (2-D matrix).
%          - Vector of X positions. If only two input arguments then this
%            can be a two column matrix of [X,Y] positions.
%          - Vector of Y positions.
%          - Filter or convolution kernel image.
%            If empty, then do not filter the image. Default is empty.
%          - Filtering function - options are:
%            'filter2' | 'filter2_fft' | 'filter2_fftc' | 'conv2' |
%            'conv_fft2'.
%            Default is 'filter2_fft'.
% Output : - Vector of pixel values at requested positions.
%          - The filtered (or unfiltered) image.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Val,Image]=ImUtil.Im.image_posval(Image,[X Y]);
%          [Val,Image]=ImUtil.Im.image_posval(Image,X,Y);
%          [Val,Image]=ImUtil.Im.image_posval(Image,X,Y,Filter,'filter2_fftc');
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<3 && size(X,2)==2)
    Y = X(:,2);
    X = X(:,1);
end

Def.Filter  = [];
Def.FiltFun = 'filter2_fft';
if (nargin<4)
    Filter  = Def.Filter;
    FiltFun = Def.FiltFun;
elseif (nargin==4)
    FiltFun = Def.FiltFun;
else
    % do nothing
end
    
if (~isempty(Filter))
    switch lower(FiltFun)
        case 'filter2'
            Image = filter2(Filter,Image,'same');
        case 'filter2_fft'
            Image = filter2_fft(Image,Filter);
        case 'filter2_fftc'
            Image = filter2_fftc(Image,Filter);
        case 'conv2'
            Image = conv2(Image,Filter,'same');
        case 'conv_fft2'
            Image = Util.external.conv_fft2(Image,Filter,'same');
        otherwise
            error('Unknown FiltFun option');
    end
end

%Ind = sub2ind(size(Image),round(Y),round(X));
% A fast replacement for sub2ind (in 2D):
Size = size(Image);
Ind  = round(Y) + (round(X)-1).*Size(1);
Val  = Image(Ind);