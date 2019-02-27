function []=measure_seeing(Image,varargin)
% Measure seeing over field of view of an image
% Package: ImUtil.pipe
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

ImageField = SIM.ImageField;
CatField   = AstCat.CatField;


DefV.ImType               = 'auto';  % 'fits' | 'hdf5' | 'mat' | 'sim' | 'auto' 
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


% extract stars
ColCell            = {'XWIN_IMAGE','YWIN_IMAGE','XPEAK_IMAGE','YPEAK_IMAGE',...
                           'X2WIN_IMAGE','Y2WIN_IMAGE','XYWIN_IMAGE',...
                           'THETA','ELONGATION',...
                           'PEAKF_VAL','PEAKF_VALTOT',...
                           'SN','SN_UNF','SN_ADD',...
                           'BACK','BACK_STD',...
                           'PEAK_GRADBACK',...
                           'FLAGS','NEAREST_SRCDIST'};
tic;
S = mextractor(S,'ColCell',ColCell,'Thresh',InPar.MinSN);
toc

SizeXY = fliplr(size(S.(ImageField)));
[ListEdge,ListCenter]=ImUtil.Im.image_blocks(SizeXY,InPar.BlockSize,10,'simple');
Nbl = size(ListEdge,1);
for Ibl=1:1:Nbl
    Flag = S.(CatField)(:,S.Col.XWIN_IMAGE)>ListEdge(Ibl,1) & ...
           S.(CatField)(:,S.Col.XWIN_IMAGE)<ListEdge(Ibl,2) & ...
           S.(CatField)(:,S.Col.YWIN_IMAGE)>ListEdge(Ibl,3) & ...
           S.(CatField)(:,S.Col.YWIN_IMAGE)<ListEdge(Ibl,4);
    X2 = S.(CatField)(Flag,S.Col.X2WIN_IMAGE);
    Y2 = S.(CatField)(Flag,S.Col.Y2WIN_IMAGE);
    XY = S.(CatField)(Flag,S.Col.XYWIN_IMAGE);
    SN = S.(CatField)(Flag,S.Col.SN);
   
    
    
end
