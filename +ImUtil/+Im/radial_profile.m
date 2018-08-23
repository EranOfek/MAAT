function [Radius,Prof,ProfErr,Npix]=radial_profile(Image,R,BinSize,MatX,MatY,ProfFun,ErrFun,varargin)
% Calculate radial profile around a point in an image.
% Package: ImUtil
% Description: Calculate radial profile around a point in an image.
% Input  : - A 2D image.
%          - Maximum radius for which to calculate profile.
%            If empty, then set to floor(ImageSize.*0.5).
%            Default is empty.
%          - Radial profile bin size. Default is 1 pix.
%          - Either a scalar of X coordinate around which to calculate the
%            radial profile or a matrix of X coordinates.
%            If empty then use image center. Default is empty.
%          - Either a scalar of Y coordinate around which to calculate the
%            radial profile or a matrix of Y coordinates.
%            If empty then use image center. Default is empty.
%          - Function handle for profile function.
%            Default is @nanmedian.
%          - Function handle for profile error.
%            Default is @nanstd.
% Output : - Radius
%          - Profile.
%          - Profile error.
%          - Number of pixels in each bin.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Radius,Prof,ProfErr,Npix]=ImUtil.Im.radial_profile(rand(10,10));
% Reliable: 2
%--------------------------------------------------------------------------

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Def.R           = [];
Def.BinSize     = 1;
Def.MatX        = [];
Def.MatY        = [];
Def.ProfFun     = @nanmedian;
Def.ErrFun      = @nanstd;
if (nargin==1)
    R           = Def.R;
    BinSize     = Def.BinSize;
    MatX        = Def.MatX;
    MatY        = Def.MatY;
    ProfFun     = Def.ProfFun;
    ErrFun      = Def.ErrFun;  
elseif (nargin==2)
    BinSize     = Def.BinSize;
    MatX        = Def.MatX;
    MatY        = Def.MatY;
    ProfFun     = Def.ProfFun;
    ErrFun      = Def.ErrFun;  
elseif (nargin==3)
    MatX        = Def.MatX;
    MatY        = Def.MatY;
    ProfFun     = Def.ProfFun;
    ErrFun      = Def.ErrFun;  
elseif (nargin==4)
    MatY        = Def.MatY;
    ProfFun     = Def.ProfFun;
    ErrFun      = Def.ErrFun;  
elseif (nargin==5)
    ProfFun     = Def.ProfFun;
    ErrFun      = Def.ErrFun;  
elseif (nargin==6)
    ErrFun      = Def.ErrFun;  
else
    % do nothing
end

[NY,NX] = size(Image);
if (isempty(R))
    R = floor(min(NX,NY).*0.5);
end

if (isempty(MatX))
    MatX = floor(NX.*0.5);
end
if (isempty(MatY))
    MatY = floor(NY.*0.5);
end

if (numel(MatX)==1)
    CenterX = MatX;
    CenterY = MatY;
    [MatX,MatY] = meshgrid((1:1:NX)-CenterX,(1:1:NY)-CenterY);
else
    [~,IX] = Util.stat.minnd(MatX);
    [~,IY] = Util.stat.minnd(MatY);
    CenterX = IX(2);
    CenterY = IY(1);
end

% BinSize = 1;
% ProfFun = @median;
% ErrFun  = @std;

X1 = max(1, floor(CenterX - R));
X2 = min(NX, ceil(CenterX + R));
Y1 = max(1, floor(CenterY - R));
Y2 = min(NY, ceil(CenterY + R));

% trim image for faster execuation
SubImage = Image(Y1:Y2,X1:X2);
SubMatX  = MatX(Y1:Y2,X1:X2);
SubMatY  = MatY(Y1:Y2,X1:X2);
SubMatR  = (SubMatX ).^2 + (SubMatY ).^2;

VecRadiusIn  = (0:BinSize:R-BinSize);
VecRadiusOut = (BinSize:BinSize:R);
% each column represent a radius bin
FlagMat       = SubMatR(:)>=VecRadiusIn.^2 & SubMatR(:)<VecRadiusOut.^2;
FlagMat       = single(FlagMat);
FlagMat(~FlagMat) = NaN;

SubImageRad   = SubImage(:).*FlagMat;
Prof    = ProfFun(SubImageRad);
ProfErr = ErrFun(SubImageRad);
Radius  = 0.5.*(VecRadiusIn+VecRadiusOut);
Npix    = nansum(FlagMat);




    