function [TileList,TileArea]=tile_the_sky(N_RA,N_Dec)
% Tile the celestial sphere
% Package: celestial.coo
% Description: Tiling the celestial sphere with approximately equal
%              area tiles.
% Input  : - Number of tiles along the celestial equator.
%          - Number of tiles along a meridian.
% Output : - List of tiles:
%            [CenterRA,CenterDec,MinRA,MaxRA,MinDec,MaxDec], in radians.
%          - Area of tiles [sr]
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [TileList,TileArea]=celestial.coo.tile_the_sky(360.*4,180.*4);
% Reliable: 2
%--------------------------------------------------------------------------

RAD     = 180./pi;
DegRA   = 360;
DegDec  = 180;

RadRA   = DegRA./RAD;
RadDec  = DegDec./RAD;

StepDec = RadDec./N_Dec;  % [radians]

VecDec  = (-0.5.*pi+0.5.*StepDec:StepDec:+0.5.*pi-0.5.*StepDec).';
N_Dec   = length(VecDec);


TileList = zeros(0,6);
TileInd  = 0;
for Idec=1:1:N_Dec
   CurrN_RA  = ceil(N_RA.*cos(VecDec(Idec) -sign(VecDec(Idec)).*0.5.*StepDec));
   StepRA    =  RadRA./CurrN_RA;    % [radians]

   VecRA     = (0+0.5.*StepRA:StepRA:2.*pi-0.5.*StepRA).';
   CurrN_RA      = length(VecRA);

   for Ira=1:1:CurrN_RA
      TileInd = TileInd + 1;
      TileList(TileInd,:) = [VecRA(Ira), VecDec(Idec), [VecRA(Ira)+0.5.*[-StepRA,+StepRA]], [VecDec(Idec)+0.5.*[-StepDec,+StepDec]]];
   end
end



TileArea = abs(celestial.coo.cel_annulus_area(TileList(:,5),TileList(:,6),'sr'))./N_RA;

