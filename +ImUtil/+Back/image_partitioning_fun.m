function [Val,X,Y,SubFun,Sub]=image_partitioning_fun(Image,BlockSize,Fun,varargin)
% Operate a function on an images in a structure array.
% Package: ImUtil.Back
% Description: Operate a function on partioned image. If a single image is
%              provided then partion the image using
%              ImUtil.Back.image_partitioning and generate a structure
%              array in which a sub image is contained in the Im field.
%              Alternativel, the input is the structure array.
%              Next operate a function on each one of the sub images and
%              return the function output.
% Input  : - An image, or a structure array of sub images (i.e., the output
%            of ImUtil.Back.image_partitioning).
%          - BlockSize [X,Y], or [X] (will be copied as [X, X]).
%            Default is empty.
%          - Function handle to operate. Default is @nanmean.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'FunAddPar' - Cell array of additional parameters to pass
%                        to the function. Default is {}.
%            'Buffer'    - Overlap between sub images. Default is 10.
% Output : - A structure array (one element per output value from Fun) of
%            matrices of function values at the sub images locations.
%          - Matrix of X coordinates of sub images centers.
%          - Matrix of Y coordinates of sub images centers.
%          - A structure array of measured values.
%          - A structure array of sub images.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Val,X,Y,SubFun,Sub]=ImUtil.Back.image_partitioning_fun(Image,[256 256],@nanmedian,'FunAddPar',{'all'});
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin<3)
    Fun = @nanmean;
end

DefV.FunAddPar            = {};
DefV.Buffer               = 10;
DefV.SecondOutPar         = false;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);




if isnumeric(Image)
    % assume Image wasn't yet partioned by  ImUtil.Back.image_partitioning
    
    if numel(BlockSize)==1
        BlockSize = [BlockSize, BlockSize];
    end
    
    [Sub,ListEdge,ListCenter] = ImUtil.Back.image_partitioning(Image,BlockSize);
    
elseif isstruct(Image)
    % do nothing
    Sub = Image;
else
    error('Unknown Image format - should be image or structure');
end
    

Nsub = numel(Sub);
SubFun = Util.struct.struct_def({'Res1','Res2'},size(Sub));
for Isub=1:1:Nsub
    SubFun(Isub).Res = Fun(Sub(Isub).Im,InPar.FunAddPar{:});
end
Npar = numel(SubFun(Isub).Res);


ListCenter = [[Sub.CenterX].', [Sub.CenterY].'];    
Nx = numel(unique(ListCenter(:,1)));
Ny = size(ListCenter,1)./Nx;

Vec = [SubFun.Res];
for Ipar=1:1:Npar
    Val(Ipar).Val = reshape(Vec(Ipar:Npar:end),Ny,Nx);
end

X   = reshape([Sub.CenterX],Ny,Nx);
Y   = reshape([Sub.CenterY],Ny,Nx);
