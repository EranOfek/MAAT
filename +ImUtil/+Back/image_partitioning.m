function [Sub,ListEdge,ListCenter]=image_partitioning(Image,BlockSize,varargin)
% Partition a 2D image into sub images.
% Package: ImUtil.Back
% Description: Partition a 2D image into sub images.
% Input  : - A single 2D image (matrix).
%          - BlockSize [X,Y], or [X] (will be copied as [X, X]).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Buffer' - Overlapping buffer. Default is 10 pix.
% Output : - A structure array in which each element contains a
%            sub image (stored in the .Im field).
%          - Four column matrix of list of blocks in image
%            [Xmin, Xmax, Ymin, Ymax].
%          - Four column matrix of list of block centers [Xcen, Ycen].
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% [Sub,ListEdge,ListCenter]=ImUtil.Back.image_partitioning(Image,[256 256]);
% Reliable: 
%--------------------------------------------------------------------------



DefV.Buffer               = 10; 
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if numel(BlockSize)==1
    BlockSize = [BlockSize, BlockSize];
end

ImageSize = fliplr(size(Image));  % [X,Y]

[ListEdge,ListCenter] = ImUtil.Im.image_blocks(ImageSize,BlockSize,InPar.Buffer,'simple');
Xmin = ListEdge(:,1);
Xmax = ListEdge(:,2);
Ymin = ListEdge(:,3);
Ymax = ListEdge(:,4);

Nblock = size(ListEdge,1);

Sub = Util.struct.struct_def({'Im','CenterX','CenterY'},Nblock,1);

for Iblock=1:1:Nblock
    Sub(Iblock).Im      = Image(Ymin(Iblock):Ymax(Iblock),Xmin(Iblock):Xmax(Iblock));
    Sub(Iblock).CenterX = ListCenter(Iblock,1);
    Sub(Iblock).CenterY = ListCenter(Iblock,2);
end

Nx = numel(unique(ListCenter(:,1)));
Ny = size(ListCenter,1)./Nx;

Sub = reshape(Sub,Ny,Nx);



  