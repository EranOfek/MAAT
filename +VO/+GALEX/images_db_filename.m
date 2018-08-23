function [FileName,FilePathNUV,FilePathFUV]=images_db_filename(DR)
% Get the database of all the GALEX images file names
% Package: VO.GALEX
% Description: Get the GALEX images file name database
% Input  : - Data release. Default is 'GR7'.
% Output : - HDF5 file name containing the catalog of GALEX
%            images. Use AstCat.loadh2astcat to read into an
%            AstCat object.
%          - mat file containing the path for the NUV images.
%          - mat file containing the path for the FUV images.
% Example: F=VO.GALEX.images_db_filename;
% Reliable: 2


if (nargin==0)
    DR = 'GR7';
end

switch lower(DR)
    case 'gr7'
        FileName = 'GALEX_GR7_Images.hdf5';
        FilePathNUV = 'GALEX_GR7_ImagesPathNUV.mat';
        FilePathFUV = 'GALEX_GR7_ImagesPathFUV.mat';
    otherwise
        error('Unknown GALEX DR option');
end

