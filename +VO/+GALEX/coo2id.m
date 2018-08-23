function [CatID,ID,Ind]=coo2id(RA,Dec,DR)
% Convert RA/Dec to GALEX image identifiers
% Package: VO.GALEX
% Description: Convert RA/Dec to GALEX image IDs.
% Input  : - Vector of J2000.0 R.A. [rad, [H M S], or sexagesimal].
%          - Vector of J2000.0 Dec. [rad, [Sign D M S], or sexagesimal].
%          - Data release. Default is 'GR7'.
% Output : - AstCat object containing the GALEX images that
%            cover the requested coordinates.
%            Each AstCat element correspond to a requested
%            coordinates.
%          - Cell array of matrices. Cell element per requested
%            coordinate. Each cell contains the GALEX image ID
%            matrix with the following columns:
%            [Vsn, Tilenum, Type, Ow, Prod, Img, Try].
%          - Cell array of index of images in the GALEX image
%            catalog and the galex image path catalog.
% Example: [CatID,ID,Ind]=VO.GALEX.coo2id('01:10:10.1','+00:20:30')
%          [CatID,ID,Ind]=VO.GALEX.coo2id([1;2],[0;1])
% Reliable: 2

RAD = 180./pi;

if (nargin==2)
    DR = 'GR7';
end

FOV_Radius = VO.GALEX.fov./RAD;   % radians


[ImagesFileName] = VO.GALEX.images_db_filename(DR);
GALEX_Images     = AstCat.loadh2astcat(ImagesFileName);  % sorted by dec
%FilePathNUV    = load2(FileNUV);
%FilePathFUV    = load2(FileFUV);
GALEX_Images    = cats.GALEX.GALEX_GR7_Images;



RA   = celestial.coo.convertdms(RA,'gH','r');
Dec  = celestial.coo.convertdms(Dec,'gD','R');
Ncoo = numel(RA);

[Res] = search_cat(GALEX_Images,RA,Dec,'SearchRad',FOV_Radius);
CatID = AstCat(Ncoo,1);
ID    = cell(Ncoo,1);
Ind   = cell(Ncoo,1);
%PathNUV = cell(Ncoo,1);
%PathFUV = cell(Ncoo,1);
for Icoo=1:1:Ncoo
    CatID(Icoo) = row_select(GALEX_Images,Res(Icoo).IndCat);
    ID{Icoo}    = col_get(CatID(Icoo),{'vsn','tileNum','Type','ow','Prod','Img','Try'});
    Ind{Icoo}   = [Res(Icoo).IndCat].';
    %PathNUV{Icoo} = FilePathNUV(Res(Icoo).IndCat);
    %PathFUV{Icoo} = FilePathFUV(Res(Icoo).IndCat);
end
