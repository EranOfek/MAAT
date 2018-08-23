function [Links,FitsName,Images,Headers]=wget_corrim(ID,Save,Filters,DataRelease)
% Get WISE coadded Atlas image from coadd_id
% Package: VO.WISE
% Description: Given a WISE field ID coadd_id for Atlas Images, get the
%              link to, and the WISE corrected image in FITS format.
%              Furthermore, read the corrected image into matlab matrix.
%              The program first check if the FITS image is exist in
%              current directory and if so, it reads the image from the disk.
%              Note that if nargout>1 then the fits file is retrieved.
% Input  : - Matrix of images ID coadd_id for Atlas Images
%          - Save FITS images to disk {'gzip' | 'y' | 'n'}, default is 'y'.
%          - Filters to retrieve, default is '1234' (all Filters).
%          - Data release {'prelim_postcryo','prelim','cryo_3band','allsky'}, default is 'allsky'.
%            If empty use default.
% Output : - Cell array of links to each fits image.
%            Rows for each field, columns for each band.
%          - Cell array of fits images name.
%          - Cell array of images matrix.
%            Rows for each field, columns for each band.
%            If cell is NaN - image can't retrieved.
%          - Cell array of images header.
%            Rows for each field, columns for each band.
%            Note that keywords are stored in: Header{:,:}.PrimaryData.Keywords
% Tested : Matlab 7.0
%     By : Ilia Labzovsky                  Oct 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Link,FN,Image,Header]=VO.WISE.wget_corrim({'0000m016_ab41'},'y','4','allsky');
%-----------------------------------------------------------------------------
RAD = 180./pi;

WISE_Server = 'http://irsa.ipac.caltech.edu/ibe';

DefSave        = 'y';
DefFilters     = '1234';
DefDataRelease = 'allsky';
if (nargin==1)
   Save         = DefSave;
   Filters      = DefFilters;
   DataRelease  = DefDataRelease;
elseif (nargin==2)
   Filters      = DefFilters;
   DataRelease  = DefDataRelease;
elseif (nargin==3)
   DataRelease  = DefDataRelease;
elseif (nargin==4)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(DataRelease))
   DataRelease = DefDataRelease;
end

switch DataRelease
case 'prelim_postcryo'
	error(fprintf('No Atlas images for prelim_postcryo \n'));
case 'prelim'
	table_name	= 'p3am_cdd';
case 'cryo_3band'
	table_name	= 'p3am_cdd';
case 'allsky'
	table_name	= '4band_p3am_cdd';
otherwise
	error(fprintf('DataRelease type not supported \n'));
end

Nband		= length(Filters);
Nim			= size(ID,1);
Links		= cell(Nim,Nband);
FitsName	= cell(Nim,Nband);
Images		= cell(Nim,Nband);
Headers		= cell(Nim,Nband);

for Iim=1:1:Nim
	% Read object ID
	coadd_id	= char(ID(Iim,1));
	coaddgrp	= coadd_id(1:2);
	coadd_ra	= coadd_id(1:4);
	URL = sprintf('%s/data/wise/%s/%s/%s/%s/%s',WISE_Server,DataRelease,table_name,coaddgrp,coadd_ra,coadd_id);
	clear table_name coaddgrp coadd_ra;
	
	for Iband=1:1:Nband
		Retrieve{Iim,Iband} = 1;
		
		ImageName = sprintf('%s-w%s-int-3.fits.gz',coadd_id,Filters(Iband));
		Links{Iim,Iband} = sprintf('%s/%s',URL,ImageName);
	
		if (nargout>1)
			switch Save
			case 'gzip'
				FitsName{Iim,Iband}= ImageName;
			otherwise
				FitsName{Iim,Iband}= ImageName(1:end-3);
			end

			%--- Check if file is already on disk ---
			if (exist(ImageName,'file')~=0 || exist(ImageName(1:end-3),'file')~=0)
				% File is already on local disk
				if (exist(ImageName,'file')==0)
					FileIsGZIP = 0;
				else
					FileIsGZIP = 1;
				end
			else
				% file is not on disk - get file
				system(sprintf('wget -q -nc %s',Links{Iim,Iband}));
				%    (sprintf('wget -q -nc %s',Links{Iim,Iband}))

				if (exist(ImageName,'file')==0)
				% Can't retriev image
					Retrieve{Iim,Iband} = 0;
				else
					Retrieve{Iim,Iband} = 1;
				end
				FileIsGZIP = 1;
			end

			switch Retrieve{Iim,Iband}
			case 1
				if (nargout>2)
					%--- read image ---
					switch FileIsGZIP
					case 1
						system(sprintf('gzip -d %s',ImageName));
					case 0 
						% do nothing
					otherwise
						error('Unknown FileIsGZIP Option');
					end 
					FitsMat = fitsread(ImageName(1:end-3));
	
					Images{Iim,Iband} = FitsMat;
					if (nargout>3)
						FitsHeader         = fitsinfo(ImageName(1:end-3));
						Headers{Iim,Iband} = FitsHeader;
					end	
				end
			 
				switch Save
				case 'gzip'
					% gzip file
					system(sprintf('gzip %s',ImageName(1:end-3)));
				case 'y'
					% do nothing - file is already ungziped
				case 'n'
					% delete FITS image (not gzip any more)
					%system(sprintf('rm %s',ImageName(1:end-3)));
					delete(sprintf('%s',ImageName(1:end-3)));
				otherwise
					error('Unknown Save Option');
				end
			case 0
				Headers{Iim,Iband} = NaN;
				Images{Iim,Iband}  = NaN;
			otherwise
				error('Unknown Retrieve Option');
			end
		end
		clear ImageName;
	end
end