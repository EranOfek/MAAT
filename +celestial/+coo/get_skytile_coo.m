function ListAllTiles=get_skytile_coo(Coo,HalfSize,N_RA,N_Dec,SubTile)
%-----------------------------------------------------------------------
% get_skytile_coo function                                    Catalogue
% Description: Assuming some sky tileing (see tile_the_sky.m) and
%              optional sub tileing for each tile, search for all the
%              tiles which their centers found within some distance
%              from a given celestial coordinates.
% Input  : - Coordinate to search [RA, Dec], in radians.
%          - Search radius in radians.
%          - Number of tiles in RA direction.
%          - Number of tiles in Dec direction.
%          - Number of subtiles in RA and Dec, within each tile.
% Output : - Tiles centers [RA, Dec], in radians.
% Tested : Matlab 7.0
%     By : Eran O. Ofek        September 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See also: tile_the_sky.m
% Example: List=get_skytile_coo([1 1.4],0.2./RAD,360,180,[5 5]);
%-----------------------------------------------------------------------

RAD       = 180./pi;
Col_RA      = 1;
Col_Dec     = 2;
Col_StepRA  = 3;
Col_StepDec = 4;
DegRA     = 360;
DegDec    = 180;

RadRA       = DegRA./RAD;
RadDec      = DegDec./RAD;

StepDec     = RadDec./N_Dec;  % [radians]

RegionN_Dec = ceil(HalfSize./StepDec);
TileDec     = (floor(Coo(2)./StepDec)+0.5).*StepDec;  % Tile center Dec

VecDec      = [TileDec-RegionN_Dec.*StepDec:StepDec:TileDec+RegionN_Dec.*StepDec].';
% remove larger than 90
I = find(abs(VecDec)<pi./2);
VecDec = VecDec(I);

ListTiles = zeros(0,4);
for Id=1:1:length(VecDec)
   CurrDec    = VecDec(Id);

                           % the upper part of the declination zone
   %CurrN_RA   = ceil(N_RA.*cos(abs(CurrDec)+StepDec.*0.5));
   CurrN_RA   = ceil(N_RA.*cos(CurrDec));
   StepRA     = RadRA./CurrN_RA;
   RegionN_RA = ceil(HalfSize./StepRA);
   TileRA     = (floor(Coo(1)./StepRA)+0.5).*StepRA;    % Tile center RA
   VecRA      = [TileRA-RegionN_RA.*StepRA:StepRA:TileRA+RegionN_RA.*StepRA].';
   % remove repeating fields
   VecRA      = unique(VecRA);

   ListTiles  = [ListTiles; [VecRA, ones(size(VecRA))*[CurrDec, StepRA, StepDec]]];


end


%---------------------------
%--- Divide to Sub tiles ---
%---------------------------
if (SubTile(1)>1 || SubTile(2)>1)
   ListAllTiles = zeros(0,2);
   N_ListTiles = size(ListTiles,1);

      
   RelSubStepRA  = 1./SubTile(1);
   RelSubStepDec = 1./SubTile(2);


   for I=1:1:N_ListTiles

      ListSubTilesRA  = [ListTiles(I,Col_RA) - 0.5.*ListTiles(I,Col_StepRA) + 0.5.*RelSubStepRA.*ListTiles(I,Col_StepRA) : RelSubStepRA.*ListTiles(I,Col_StepRA) : ListTiles(I,Col_RA) + 0.5.*ListTiles(I,Col_StepRA) - 0.5.*RelSubStepRA.*ListTiles(I,Col_StepRA)];

      ListSubTilesDec = [ListTiles(I,Col_Dec) - 0.5.*ListTiles(I,Col_StepDec) + 0.5.*RelSubStepDec.*ListTiles(I,Col_StepDec) : RelSubStepDec.*ListTiles(I,Col_StepDec) : ListTiles(I,Col_Dec) + 0.5.*ListTiles(I,Col_StepDec) - 0.5.*RelSubStepDec.*ListTiles(I,Col_StepDec)];

      [X,Y] = meshgrid(ListSubTilesRA,ListSubTilesDec);

      ListAllTiles = [ListAllTiles; [Util.array.mat2vec(X), Util.array.mat2vec(Y)]];
   end
else
   ListAllTiles = ListTiles;
end
   

%----------------------------------
%--- Convert RA to 0..2pi range ---
%----------------------------------
I = find(ListAllTiles(:,Col_RA)<0);
ListAllTiles(I,Col_RA) = 2.*pi + ListAllTiles(I,Col_RA);


Dist = sphere_dist(Coo(Col_RA),Coo(Col_Dec),ListAllTiles(:,Col_RA),ListAllTiles(:,Col_Dec));

CorrectedHalfSize = (HalfSize + max([N_Dec.*2./360;N_RA./360]).*2./min(SubTile)./RAD./sqrt(2) );


I = find(Dist<=CorrectedHalfSize);

ListAllTiles = ListAllTiles(I,:);


