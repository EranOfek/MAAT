function [SectorBack,SecSub,FiltSecSub]=sectorial_median_filter(Image,X,Y,varargin)
% Apply a sectorial median filter to an image.
% Package: ImUtil.Im
% Description: Given an image an some coordinates, calculate he sectorial
%              median filter. For each point, a sector on this point,
%              relative to the central coordinates is calculated and the
%              median of this sector is replaing the position value.
% Input  : - An image.
%          - X coordinate of image center relative to which all sectors
%            are calculated.
%          - Y coordinate of image center.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'DR'    - Sector semi width. Default is 1 pix.
%            'Dphi'  - Sector semi opening angle. Default is 30 deg.
%            'ExRad' - Radius around pixel from which to exclude from the
%                      sector median calculation. Default is 2 pix.
%            'Filter'- A matrix of filter to apply to the image after
%                      sectorial median is subtracted.
%                      Default is Kernel2.gauss(1.5,1.5).
% Output : - Image of the sectorial median in each point, relative to the
%            central coordinates.
%          - The original image from which the sectorial medain is
%            subtracted.
%          - The median sectorial subtracted image filtered with the input
%            filter.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.DR                  = 1;
DefV.Dphi                = 30;
DefV.ExRad               = 2;
DefV.Filter              = Kernel2.gauss(1.5,1.5);
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Dphi = InPar.Dphi./RAD;

SizeIm = size(Image);
[MatX, MatY] = meshgrid((1:1:SizeIm(2))-X,(1:1:SizeIm(1))-Y);
% distance from image center
MatR = sqrt(MatX.^2 + MatY.^2);
% angle relative to image center
MatPhi = atan2(MatY,MatX);
% convert to 0..2pi range
MatPhi(MatPhi<0) = 2.*pi + MatPhi(MatPhi<0);

% for each Xp, Yp
SectorBack = zeros(SizeIm(1),SizeIm(2));
for Xp=1:1:SizeIm(2)
    for Yp=1:1:SizeIm(1)
        
        [MatXp,MatYp] = meshgrid((1:1:SizeIm(2))-Xp,(1:1:SizeIm(1))-Yp);
        MatRp         = sqrt(MatXp.^2 + MatYp.^2);
        DX  = Xp - X;
        DY  = Yp - Y;
        R   = sqrt(DX.^2 + DY.^2);
        Phi = atan2(DY,DX);
        Phi(Phi<0) = 2.*pi + Phi(Phi<0);
        RelPhi = MatPhi - Phi;

        Ipix = find(MatR>=(R-InPar.DR) & MatR<=(R+InPar.DR) & ...
                    (abs(RelPhi)>=(2.*pi-Dphi) | abs(RelPhi)<Dphi) & ...
                    MatRp>InPar.ExRad);
         SectorBack(Yp,Xp) = median(Image(Ipix));
    end
end

if (nargout>1)
    % subtract sectorial median
    SecSub = Image - SectorBack;

    if (nargout>2)
        % filter
        FiltSecSub = filter2(InPar.Filter,SecSub);
    end
end


     


