function [D]=wmedian_im(Images,Err,Perct)
%--------------------------------------------------------------------------
% wmedian_im function                                              General
% Description: Weighted median on a set of images in a cube, where the
%              image index is in the 1st dimension.
% Input  : - Cube of images in which the image index is the 1st dimension.
%          - Vector of error, per image.
%          - Percentile. Defaiult is 0.5 (i.e., median).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - Weighted median image. Note that sometimes the weighted
%            median is undefined.
% License: GNU general public license version 3
% Tested : Matlab R2015a
%     By : Eran O. Ofek                    Aug 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [D]=wmedian_im(rand(5,10,10),rand(5,1))
% Reliable: 
%--------------------------------------------------------------------------

why median is undefined?


if (nargin==2),
    Perct = 0.5;
end

% transform the cube into a matrix [ImageInd X PixelInd]
Size   = size(Images);
Nim    = Size(1);
Npix   = Size(2).*Size(3);
Images = reshape(Images,[Nim Npix]);

[SV,SI] = sort(Images);
CSW     = cumsum(1./Err(SI).^2);  % cumsum of the weights
CSWnorm = bsxfun(@times,CSW,1./CSW(end,:));
D       = zeros(Npix,1);
for Ipix=1:1:Npix,
    D(Ipix) = interp1(CSWnorm(:,Ipix),SV(:,Ipix),Perct,'linear');
end
D = reshape(D,Size(2:3));





