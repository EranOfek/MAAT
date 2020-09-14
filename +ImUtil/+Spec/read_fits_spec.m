function Spec=read_fits_spec(varargin)
%--------------------------------------------------------------------------
% read_fits_spec function                                           ImSpec
% Description: Get a spectrum from a FITS file.
%              Currently support on CTYPE1='linear'.
% Input  : * see fitsread for options
%            e.g., FileName,Extension,...
%          - Spectrum [Wavelength [A], Flux, Pixel]
% Tested : Matlab 7.8
%     By : Eran O. Ofek                   Dec 2009
%    URL : http://weizmann.ac.il/home/eran/matlab/
% Example: Spec=read_fits_spec('ptf09uj.all.fits');
% see also: read_fits_spec.m
% Reliable: 2
%--------------------------------------------------------------------------

Spec = fitsread(varargin{:}).';
N    = length(Spec);
VecPix = (1:1:N).';


KeywordVal=get_fits_keyword(varargin{1},{'CRVAL1','CDELT1','CRPIX1','CTYPE1','CD1_1'});
CRVAL1 = KeywordVal{1};
CDELT1 = KeywordVal{2};
CRPIX1 = KeywordVal{3};
CTYPE1 = KeywordVal{4};
CD1_1  = KeywordVal{5};

if (isnan(CDELT1)==1)
   CDELT1 = CD1_1;
end


switch lower(deblank(CTYPE1))
 case 'linear'
    Wave = CRVAL1 + (VecPix - CRPIX1).*CDELT1;

    Spec = [Wave, Spec, VecPix];
 
 otherwise
    error('Unsported CTYPE1 value');
end
