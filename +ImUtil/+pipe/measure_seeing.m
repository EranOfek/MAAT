function [Res]=measure_seeing(Image,varargin)
% Measure seeing over field of view of an image
% Package: ImUtil.pipe
% Description: Measure seeing in an image
% Input  : - Image, in FITS, HDF5, mat or SIM format.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ImType'   - One of the following:
%                         'fits', 'hdf5', 'mat', 'sim', 'auto'.
%                         Default is 'auto' (automatic identification).
%            'Section'  - Central part of image for which to calc seeing.
%                         1 - for all image, <1 for central part.
%                         Default is 1.
%            'MinSN'    - FFU
%            'MaxSN'    - FFU
%            'BlockSize'- FFU
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=ImUtil.pipe.measure_seeing(S(1),'Section',0.1)
% Reliable: 2
%--------------------------------------------------------------------------

ImageField = SIM.ImageField;
CatField   = AstCat.CatField;


DefV.ImType               = 'auto';  % 'fits' | 'hdf5' | 'mat' | 'sim' | 'auto' 
DefV.Section              = 1; % or <1 number
DefV.MinSN                = 10;
DefV.MaxSN                = 1000;
DefV.BlockSize            = [256 256];
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


switch lower(InPar.ImType)
    case 'auto'
        if (SIM.issim(Image))
            InPar.ImType = 'sim';
        else
            if (isnumeric(Image))
                InPar.ImType = 'mat';
            else
                if (ischar(Image))
                    if strcmp(Image(end-5:end),'.fits')
                        InPar.ImType = 'fits';
                    else
                        InPar.ImType = 'hdf5';
                    end
                end
            end
        end
end
                        

% read image
switch lower(InPar.ImType)
    case 'fits'
        S = FITS.read2sim(Image);
    case 'hdf5'
        error('hdf5 format is not supported yet');
    case 'mat'
        S = SIM;
        S.(ImageField) = Image;
    case 'sim'
        S = Image;
    case 'auto'
        error('Can not identify image type automatically');
end

SizeXY = fliplr(size(S.(ImageField)));

if (InPar.Section==1)
    % do nothing
else
    Cor    = round(0.5.*SizeXY + 0.5.*InPar.Section.* SizeXY.*[1;-1]);
    CCDSEC = [Cor(2,1) Cor(1,1) Cor(2,2) Cor(1,2)];
    S      = trim_image(S,CCDSEC);
end
    


% extract stars
ColCell            = {'XWIN_IMAGE','YWIN_IMAGE','XPEAK_IMAGE','YPEAK_IMAGE',...
                           'X2WIN_IMAGE','Y2WIN_IMAGE','XYWIN_IMAGE',...
                           'THETA','ELONGATION',...
                           'PEAKF_VAL','PEAKF_VALTOT',...
                           'SN','SN_UNF','SN_ADD',...
                           'BACK','BACK_STD',...
                           'PEAK_GRADBACK',...
                           'FLAGS','NEAREST_SRCDIST'};

[B,E] = ImUtil.Im.background_fit(S.Im);
S.BackIm = B;
S.ErrIm  = E;
S = mextractor(S,'ColCell',ColCell,'Thresh',InPar.MinSN);

[CG,HWHM]=S.curve_growth_psf;
FWHM = 2.*HWHM;   % pix

X2 = S.(CatField)(:,S.Col.X2WIN_IMAGE);
Y2 = S.(CatField)(:,S.Col.Y2WIN_IMAGE);
XY = S.(CatField)(:,S.Col.XYWIN_IMAGE);
SN = S.(CatField)(:,S.Col.SN);

Res.PSF_FWHM_pix = FWHM;


% SizeXY = fliplr(size(S.(ImageField)));
% [ListEdge,ListCenter]=ImUtil.Im.image_blocks(SizeXY,InPar.BlockSize,10,'simple');
% Nbl = size(ListEdge,1);
% for Ibl=1:1:Nbl
%     Flag = S.(CatField)(:,S.Col.XWIN_IMAGE)>ListEdge(Ibl,1) & ...
%            S.(CatField)(:,S.Col.XWIN_IMAGE)<ListEdge(Ibl,2) & ...
%            S.(CatField)(:,S.Col.YWIN_IMAGE)>ListEdge(Ibl,3) & ...
%            S.(CatField)(:,S.Col.YWIN_IMAGE)<ListEdge(Ibl,4);
%     X2 = S.(CatField)(Flag,S.Col.X2WIN_IMAGE);
%     Y2 = S.(CatField)(Flag,S.Col.Y2WIN_IMAGE);
%     XY = S.(CatField)(Flag,S.Col.XYWIN_IMAGE);
%     SN = S.(CatField)(Flag,S.Col.SN);
%    
%     
% end
