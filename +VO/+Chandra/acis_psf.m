function [PSF,PSFxVec,PSFyVec]=acis_psf(Energy,DetX,DetY,varargin)
% Read and interpolate the Chandra ACIS-S/ACIS-I PSF.
% Package: VO.Chandra
% Description: Read and interpolate the Chandra ACIS-S/ACIS-I PSF from
%              the Chandra CalDB files.
% Instellation: 1. Install the Chandra CalDB directory including the 2d
%                  PSFs.
%               2. Modify the DefV.CalDBdir to direct to the Chandra CalDB
%                  directory.
% Input  : - Energy [keV]
%          - Detector X position (DetX) [pix].
%            The center of the detector is at 4096.5
%          - Detector Y position (DetY) [pix].
%            The center of the detector is at 4096.5
%          * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'Camera'  - {'ACIS-S','ACIS-I'}. Default is 'ACIS-S'.
%            'FileName'- PSF file name in the CalDB directory.
%                        Default is 1 for 'aciss1998-11-052dpsf1N0002.fits'.
%                        1 - in the range of -10:1:+10 arcmin
%                            pix scale 6 micron, 256x256 PSF image, 580 Mb
%                        2 - in the range of -25:5:+10 arcmin
%                            pix scale 12 micron, 512x512 PSF image, 204 Mb
%                        3 - in the range of -6:1:+6 arcmin
%                            pix scale 2 micron, 512x512 PSF image, 890 Mb
%                        4 - in the range of -1:1:+1 arcmin
%                            pix scale 1 micron, 512x512 PSF image, 46 Mb
%            'SaveFITS'- A string of a FITS file name in which to save
%                        the PSF. If empty matrix then do not save file.
%                        Default is empty matrix.
% Output : - 2D matrix containing the Chandra PSF at the requested energy
%            and position.
%          - Vector of X coordinates corresponding to the PSF image in
%            units of pixels (pixel scale is 0.492"/pix).
%          - Vector of Y coordinates corresponding to the PSF image in
%            units of pixels (pixel scale is 0.492"/pix).
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Jan 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [PSF,X,Y]=VO.Chandra.acis_psf(1,4090,4090,'FileName',1);
% Reliable: 2
%--------------------------------------------------------------------------

CalDBrel = 'data/chandra/acis/2d_psf';

DefV.InterpMethod = 'nearest';
DefV.Camera       = 'ACIS-S';
DefV.FileName     = 1;
DefV.CalDBdir     = '/home/eran/matlab/data/private/Chandra/CalDB';
DefV.SaveFITS     = [];
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

CalDBdir = sprintf('%s%s%s',InPar.CalDBdir,filesep,CalDBrel);

switch lower(InPar.Camera)
    case {'acis-s','aciss'}
        CameraFN = 'aciss';
    case {'acis-i','acisi'}
        CameraFN = 'acisi';
    otherwise
        error('Unknown Camera option');
end

if (ischar(InPar.FileName))
    % FileName is already a string
else
    % FileName is an integer
    InPar.FileName = sprintf('%s1998-11-052dpsf%1dN0002.fits',CameraFN,InPar.FileName);
end 

FullPathFileName = sprintf('%s%s%s',CalDBdir,filesep,InPar.FileName);
Im       = fitsread(FullPathFileName);
SizeIm   = size(Im);
Im2      = fitsread(FullPathFileName,'BinTable',2);
% CTYPE1P - parameter name
% CRPIX1P - Center pixel of PSF image in PSF image pixels
% CRVAL1P - Center pixel of PSF image in AXAF Focal Plane 
% CD1_1P  - AXAF focal plane pixels per PSF image pixel
% 1 - PSFX
% 2 - PSFY
% 3 - defocus
% 4 - energy
% 5 - DETX
% 6 - DETY

Keys     = {'CTYPE1P','CRPIX1P','CRVAL1P','CD1_1P',...
            'CTYPE2P','CRPIX2P','CRVAL2P','CD2_2P',...
            'CRVAL5P','CD5_5P','CRPIX5P',...
            'CRVAL6P','CD6_6P','CRPIX6P'};
[~,~,KeyS]=FITS.get_keys(FullPathFileName,Keys,HDUnum,true);
%[~,KeyS] = get_fits_keyword(FullPathFileName,Keys);


PSFxVec    = (-(KeyS.CRPIX1P-1):1:KeyS.CRPIX1P).'.*KeyS.CD1_1P;
PSFyVec    = (-(KeyS.CRPIX2P-1):1:KeyS.CRPIX2P).'.*KeyS.CD2_2P;
DefocusVec = 0;
EnergyVec  = Im2{2};
DetxVec    = ((KeyS.CRVAL5P - KeyS.CD5_5P.*(KeyS.CRPIX5P-1)):KeyS.CD5_5P:(KeyS.CRVAL5P + KeyS.CD5_5P.*(KeyS.CRPIX5P-1))).';
DetyVec    = ((KeyS.CRVAL6P - KeyS.CD6_6P.*(KeyS.CRPIX6P-1)):KeyS.CD6_6P:(KeyS.CRVAL6P + KeyS.CD6_6P.*(KeyS.CRPIX6P-1))).';

% remove the defocus dimension
Im = squeeze(Im);

PSF = interpn(PSFxVec,PSFyVec,EnergyVec,DetxVec,DetyVec,...
              Im,...
              PSFxVec,PSFyVec,Energy,DetX,DetY,...
              InPar.InterpMethod);
    

if (~isempty(InPar.SaveFITS))
    % save FITS
    HeaderInfo = FITS.cellhead_addkey([],'COMMENT','','Chandra PSF created by chandra_psf.m',...
                                         'Camera',InPar.Camera,'Chandra camera',...
                                         'Interp',InPar.InterpMethod,'PSF interpolation method',...
                                         'CalDB',InPar.CalDBdir,'CalDB dir',...
                                         'File',InPar.FileName,'Original PSF file');

    FITS.write(PSF,InPar.SaveFITS,'Header',HeaderInfo);
end
    