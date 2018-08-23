function All=find_within_radius_mat(Image,X,Y,Radius)
% Find points within a radius from a list of coordinates.
% Package: ImUtil.Im
% Description: Given a list of coordinates within an array, return for each
%              coordinate a vector (within a matrix) of indices of
%              the points in the image
%              that are within a given box with a given radius from the
%              coordinate.
%              Also return a corresponding vector of distances squared of
%              each point in the image from the coordinate.
%              Also return the X and Y coordinates of each point and the
%              value in the image.
%              See also find_within_radius_cell.m (somewhat faster).
% Input  : - An image.
%          - A vector of X coordinates around to search for the points
%            within the radius.
%          - A vector of Y coordinates around to search for the points
%            within the radius.
%          - Radius [pix].
% Output : - A structure array with the following fields:
%            .MatI - A matrix of indices of the points within the image.
%                    Each column corresponds to a requested coordinates.
%                    All the points are within a box which is width is
%                    1+2*radius centered on the requested coordinate.
%            .MatX - Like .MatI, but for the X coordinates of the points.
%            .MatY - Like .MatI, but for the Y coordinates of the points.
%            .MatR2- Like .MatI, but for the radius squared of each point.
%            .MatV - Like .MatI, but for the image value of each point.
% See also: find_within_radius_cell.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: All = ImUtil.Im.find_within_radius_mat(Image,X,Y,20);
% Reliable: 
%--------------------------------------------------------------------------

error('bug!')

Nsrc = numel(X);

Size = size(Image);
Npt  = ((Radius+1).*2).^2;
All.MatX  = nan(Npt,Nsrc);
All.MatY  = nan(Npt,Nsrc);
All.MatR2 = nan(Npt,Nsrc);
All.MatV  = nan(Npt,Nsrc);
All.MatI  = zeros(Npt,Nsrc,'uint32');
for Isrc=1:1:Nsrc,
    MinI = max(floor(Y(Isrc) - Radius),1);
    MaxI = min(ceil(Y(Isrc) + Radius), Size(1));
    
    MinJ = max(floor(X(Isrc) - Radius),1);
    MaxJ = min(ceil(X(Isrc) + Radius), Size(2));
    
    [MatX,MatY] = meshgrid((MinJ:MaxJ),(MinI:MaxI));
    %Rad2{Isrc}   = (MatX-X(Isrc)).^2 + (MatY-Y(Isrc)).^2;
    
    %[Npt, numel(MatX)]
    N = numel(MatX);
    All.MatX(1:N,Isrc)  = MatX(:);
    All.MatY(1:N,Isrc)  = MatY(:);
    All.MatR2(1:N,Isrc) = MatX(:).^2 + MatY(:).^2;
    %All.MatI(1:N,Isrc)  = sub2ind(Size,MatY(:),MatX(:));
    % A fast replacment for sub2ind (in 2D):
    All.MatI(1:N,Isrc) = MatY(:) + (MatX(:)-1).*Size(1);
    All.MatV(1:N,Isrc)  = Image(All.MatI(1:N,Isrc));
end
    
