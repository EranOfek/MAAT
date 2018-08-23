function [Section,FlagIn,CenterSection]=best_section(ImageSizeXY,CenterXY,SectionHalfSize,MinDist)
% Select best ovelapping section containing a point in a set of images
% Package: ImUtil.Im
% Description: Given images size and list of position within each image
%              which correspond to the same physical coordinate, select a
%              section in each image that contains the requested position.
%              The sections are chosen such that the requested position
%              have some minimal distance from the image boundry. Images in
%              which such section can not be found are flagged.
%              This is useful in order to find a trim region in images that
%              are completly overlapping (i.e., the requested position have
%              the same position within the section), handave the same
%              size.            
%             Note that this function is not optimal as it is not
%             maximizing the number of valid sections.
% Input  : - A two column matrix of [X, Y] image sizes.
%          - A two column matrix of [X, Y] centeres in the images.
%            The sections will be selected to contain these locations.
%          - The section half size [X,Y]. Default is [511 511].
%          - Minimum distance of position from section boundry.
%            Default is 64.
% Output : - A four column matrix of all sections (for all images)
%            [Xstart, Xend, Ystart, Yend].
%          - A logical flag indicating that the section is within image
%            boundries.
%          - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ImageSizeXY=[2024 2024; 2024 2024; 2024 2000; 2024 2024];
%          CenterXY=[100 100;10 500;1000 1000 ;1900 1900];
%          [Section,FlagIn]=ImUtil.Im.best_section(ImageSizeXY,CenterXY);
% Reliable: 2
%--------------------------------------------------------------------------

Def.SectionHalfSize = [511 511];
Def.MinDist         = 64;

if (nargin==2)
    SectionHalfSize = Def.SectionHalfSize;
    MinDist         = Def.MinDist;
elseif (nargin==3)
    MinDist         = Def.MinDist;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments: best_section(ImageSizeXY,CenterXY,[SectionHalfSize,MinDist])');
end

% Select images in which target is within image and at least MinDist from
% boundries
FlagIn = CenterXY(:,1)>MinDist & CenterXY(:,2)>MinDist & (ImageSizeXY(:,1)-CenterXY(:,1))>MinDist & (ImageSizeXY(:,2)-CenterXY(:,2))>MinDist;

MinX = min(CenterXY(FlagIn,1))-1;
MinY = min(CenterXY(FlagIn,2))-1;

CenterSection = [min(MinX,SectionHalfSize(1)+1), min(MinY,SectionHalfSize(2)+1)];

SectionHalfSizeX1 = min(MinX,SectionHalfSize(1));
SectionHalfSizeY1 = min(MinY,SectionHalfSize(2));
SectionHalfSizeX2 = SectionHalfSize(1).*2 - SectionHalfSizeX1;
SectionHalfSizeY2 = SectionHalfSize(2).*2 - SectionHalfSizeY1;


% Define section boundries
Section = [CenterXY(:,1) - SectionHalfSizeX1, CenterXY(:,1) + SectionHalfSizeX2, ...
           CenterXY(:,2) - SectionHalfSizeY1, CenterXY(:,2) + SectionHalfSizeY2];

% check that section boundries are within image size
FlagIn = FlagIn & Section(:,2)<=ImageSizeXY(:,1) & Section(:,4)<=ImageSizeXY(:,2);

