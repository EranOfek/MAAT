function [CellInd,CellRad2]=find_within_radius_cell(Size,X,Y,Radius,Circle)
% Find points within a radius from a list of coordinates.
% Package: ImUtil.Im
% Description: Given a list of coordinates within an array, return for each
%              coordinate a vector (within a cell array) of indices of
%              the points in the image
%              that are within a given radius from the coordinate.
%              Also return a corresponding vector of distances squared of
%              each point in the image from the coordinate.
%              See also find_within_radius_mat.m, but this function is
%              usually faster.
% Input  : - An array size [Y, X] (e.g., size(Image)).
%          - A vector of X coordinates around to search for the points
%            within the radius.
%          - A vector of Y coordinates around to search for the points
%            within the radius.
%          - Radius [pix].
%          - If true, then will return only points within the radius.
%            If false, then will return all the points within a square of
%            size: radius*2+1.
%            Default is false (somewhat faster).
% Output : - A cell array of vectors. Each vector contains the indices of
%            points that are found within the radius from the point.
%            Each cell element corresponds to a coordinate in the X,Y
%            vectors.
%          - A cell array of vectors. Each vector contains the distance
%            squares of the points in the image from the requested
%            coordinate.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [CI,CR2]=ImUtil.Im.find_within_radius_cell([1024 1024],X,Y,20);
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin==4)
    Circle = false;
end

Radius2 = Radius.^2;
Nsrc = numel(X);
%[MatX,MatY]=meshgrid((1:1:Size(2)),(1:1:Size(1)));
CellInd  = cell(Nsrc,1);
CellRad2 = cell(Nsrc,1);

% define template of MatX and MatY
Npt = 1+Radius.*2;
[TempMatX,TempMatY] = meshgrid((1:Npt),(1:Npt));

for Isrc=1:1:Nsrc
    % define the indices of the image stamp around the source
    % Y coordinates of image stamp
    MinI = max(floor(Y(Isrc) - Radius),1);
    MaxI = min(ceil(Y(Isrc) + Radius), Size(1));
    % X coordinates of image stamp
    MinJ = max(floor(X(Isrc) - Radius),1);
    MaxJ = min(ceil(X(Isrc) + Radius), Size(2));
    % Define the coordinate system within the stamp
    
    if (MinJ>1 && MinI>1 && MaxJ<Size(2) && MaxI<Size(1))
        % stamp is well within image
        % use previous stamp if possible
        MatX = TempMatX + MinJ - 1;
        MatY = TempMatY + MinI - 1;
    else
        [MatX,MatY] = meshgrid((MinJ:MaxJ),(MinI:MaxI));
    end
    CellRad2{Isrc}   = (MatX-X(Isrc)).^2 + (MatY-Y(Isrc)).^2;
    
    if (Circle)
        % Search in circle
        Flag      = CellRad2{Isrc}<Radius2;
        CellInd{Isrc} = sub2ind(Size,MatY(Flag),MatX(Flag));
        CellRad2{Isrc} = CellRad2{Isrc}(Flag);
    else
        % search in box
        %CellInd{Isrc} = sub2ind(Size,MatY(:),MatX(:));
        % A fast replacment for sub2ind (in 2D):
        CellInd{Isrc} = MatY(:) + (MatX(:)-1).*Size(1);
    end
end
    
