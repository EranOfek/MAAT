function [PathNUV,PathFUV]=ind2path(Ind,DR)
% Convert index of the GALEX images DB file to GALEX image path
% Package: VO.GALEX
% Description: Convert index in the GALEX images file to
%              GALEX images path.
% Input  : - Index.
%          - Data relaese. Default is 'GR7'.
% Output : - GALEX NUV images path.
%          - GALEX FUV images path.
% Reliable: 2

if (nargin==1)
    DR = 'GR7';
end

[~,FileNUV,FileFUV] = VO.GALEX.images_db_filename(DR);
FilePathNUV    = Util.IO.load2(FileNUV);
FilePathFUV    = Util.IO.load2(FileFUV);
PathNUV        = FilePathNUV(Ind);
PathFUV        = FilePathFUV(Ind);
