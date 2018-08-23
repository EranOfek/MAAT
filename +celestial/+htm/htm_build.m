function [HTM,LevList]=htm_build(Level,IsClassHTM)
% Build Hierarchical Triangular Mesh (HTM) structure
% Package: celestial.htm
% Description: Build Hierarchical Triangular Mesh (HTM) structure.
%              This structure can be use for fast searches of data
%              in catalogs on a sphere.
% Input  : - The number of levels in the HTM structure.
%          - A flag indicating if to return a ClassHTM object. Default is
%            false.
% Output : - The HTM structure array with the follwoing fields.
%            .level  - Level depth  index (0 for the first level).
%            .coo    - Coordinates of triangular mesh [Long, Lat] in
%                      radians. (3x2 matrix).
%                      The coordinates are ordered such that the
%                      right-hand rule is pointing toward the
%                      center of the polygon.
%            .cosd   - Cosine directions of triangular mesh.
%                      3x3 matrix in which each line corresponds to
%                      each vertex of the triangle.
%            .id     - Triangle id. This is a vector in which the
%                      number of elements equal to the number of levels.
%                      The first level is between 0 to 7 and all
%                      the other levels are between 0 to 3.
%            .father - The index of the father.
%                      If empty matrix then there is no father.
%            .son    - Vector of indices of all sons.
%                      If empty matrix then there are no sons.
%          - A structure array of list of levels (LevList).
%            Number of elements corresponds to the number of levels.
%            The structure contains the following fields:
%            .level - Level depth  index (0 for the first level).
%            .ptr   - Vector of indices of all elements in HTM
%                     which are in this level.
%            .side  - Length of side of triangles in level [radians].
% Tested : Matlab 7.11
%     By : Eran O. Ofek                      July 2011
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [HTM,LevList]=celestial.htm.htm_build(4);
% Reliable: 2
%------------------------------------------------------------------------------

if (nargin<2)
    IsClassHTM = false;
end

% check if data exist
% LoadFromFile = false;
% if (nargin==1)
%     DataFileName = sprintf('HTM_%d.mat',Level);
% 
%     if (exist(DataFileName,'file')>0)
%         LoadFromFile = true;
%     end
% end
% 
% if (LoadFromFile)
%     load(DataFileName);
% else
    
    
% build northern zero-level HTM
Ind = 0;

if (IsClassHTM)
    HTM = ClassHTM;
else
    HTM = Util.struct.struct_def({'level','id','coo','cosd','center_coo','center_cosd','father','son','cat'},1,1);
end
for I=1:1:4
   Ind = Ind + 1;
   HTM(Ind).level    = 0;
   HTM(Ind).coo      = [0 0; pi./2 0; 0 pi./2];
   HTM(Ind).coo(:,1) = HTM(Ind).coo(:,1) + pi./2.*(I-1);
   [CD1, CD2, CD3]   = celestial.coo.coo2cosined(HTM(Ind).coo(:,1),HTM(Ind).coo(:,2));
   HTM(Ind).cosd     = [CD1, CD2, CD3];  
   %HTM(Ind).center_cosd = mean(HTM(Ind).cosd);
   %[LonC,LatC]       = celestial.coo.cosined2coo(HTM(Ind).center_cosd(1),HTM(Ind).center_cosd(2),HTM(Ind).center_cosd(3));
   %HTM(Ind).center_coo = [LonC, LatC];
   HTM(Ind).id       = Ind-1;
   HTM(Ind).father   = [];
   %HTM(Ind).brother  = [1:1:8];
   HTM(Ind).son      = [];
   HTM(Ind).cat      = [];
   [PoleLong,PoleLat]  = celestial.htm.polysphere_poles(HTM(Ind).coo(:,1), HTM(Ind).coo(:,2));
   HTM(Ind).PolesCoo   = [PoleLong(:), PoleLat(:)];
end


% build south hemisphere zero-level HTM
for I=1:1:4
   Ind = Ind + 1;
   HTM(Ind).level    = 0;
   HTM(Ind).coo      = [0 0; 0 -pi./2; pi./2 0];
   HTM(Ind).coo(:,1) = HTM(Ind).coo(:,1) + pi./2.*(I-1);
   [CD1, CD2, CD3]   = celestial.coo.coo2cosined(HTM(Ind).coo(:,1),HTM(Ind).coo(:,2));
   HTM(Ind).cosd     = [CD1, CD2, CD3];  
   %HTM(Ind).center_cosd = mean(HTM(Ind).cosd);
   %[LonC,LatC]       = celestial.coo.cosined2coo(HTM(Ind).center_cosd(1),HTM(Ind).center_cosd(2),HTM(Ind).center_cosd(3));
   %HTM(Ind).center_coo = [LonC, LatC];
   HTM(Ind).id       = Ind-1;
   HTM(Ind).father   = [];
   %HTM(Ind).brother  = [1:1:8];
   HTM(Ind).son      = [];
   HTM(Ind).cat      = [];
   HTM(Ind).cat      = [];
   [PoleLong,PoleLat]  = celestial.htm.polysphere_poles(HTM(Ind).coo(:,1), HTM(Ind).coo(:,2));
   HTM(Ind).PolesCoo   = [PoleLong(:), PoleLat(:)];
end

LevList(1).level = 0;
LevList(1).ptr   = (1:1:8);
LevList(1).side  = pi./2;


if (Level==0)
   % stop
else
   % build HTM recursively
   [HTM,LevList] = celestial.htm.htm_build_son(HTM,LevList,Level,Ind);
end
 
%end