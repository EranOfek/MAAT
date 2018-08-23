function swarp(ImList,varargin)
%--------------------------------------------------------------------------
% swarp function                                                  ImAstrom
% Description: Run SWarp.
% Input  : - Set of images. See create_list.m for options.
%          - Alternatively, this can be a SIM or a structure array of
%            images. In this case the images will be saved as a temporary
%            FITS files.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'OutIm'    - Output coadd image name. Default is 'coadd.fits'.
%            'OutWeightIm'- Output weight coadd image name.
%                         Default is 'coadd.weight.fits'.
%            'RA'       - Center RA parameter. Radians, [H M S] or
%                         sexagesimal string. Default is [].
%            'Dec'      - Center Dec parameter. Radians, [sign D M S] or
%                         sexagesimal string. Default is [].
%            'SWarpPars'- String or cell array of pairs of ...,key,val,...
%                         SWarp parameters.
%            'SWarpProg'- SWarp program name (e.g., 'swarp','swarp-mp').
%                         Default is 'swarp'
%            'SWarpDir' - SWarp directory (absolute or relative to this
%                         function). Default is '/usr/local/bin/'.
%            'ConfigFile'- Swarp configuration file. Default is ''.
%            'Weight'   - Vector of weighst per image, that will be used
%                         to updae the FSCALE_KEYWORD header keyword.
%                         Default is equal weights.
% Output : null 
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% ImUtil.Im.swarp(ImList,'SWarpPars',{'PIXELSCALE_TYPE','MANUAL','PIXEL_SCALE','1.0','RESAMPLING_TYPE','LANCZOS2',OVERSAMPLING','3,3','COMBINE_TYPE','WEIGHTED','IMAGE_SIZE','2250,2250'});
% ImUtil.Im.swarp(ImList,'RA',1.589,'Dec',0.191,'SWarpPars',{'PIXELSCALE_TYPE','MANUAL','PIXEL_SCALE','1.0','RESAMPLING_TYPE','LANCZOS2',OVERSAMPLING','3,3','COMBINE_TYPE','WEIGHTED','IMAGE_SIZE','2250,2250'});
% % to generate background subtracted resampled images:
% ImUtil.Im.swarp(ImList,'RA',1.589,'Dec',0.191,'SWarpPars',{'COMBINE','N','PIXELSCALE_TYPE','MANUAL','PIXEL_SCALE','1.0','RESAMPLING_TYPE','LANCZOS2','OVERSAMPLING','3,3','COMBINE_TYPE','WEIGHTED','IMAGE_SIZE','2250,2250'});
% Reliable: 2
%--------------------------------------------------------------------------
import celestial.coo.*
import Util.string.*

DefV.OutIm          = [];
DefV.OutWeightIm    = [];
DefV.SWarpProg      = 'swarp';   % {'swarp'|'swarp-mp'}
DefV.SWarpDir       = '/usr/local/astro-sw/bin/'; %'/usr/local/bin/';
DefV.ConfigFile     = '';
DefV.RA             = [];
DefV.Dec            = [];
DefV.SWarpPars      = {};
DefV.Weight         = [];
DefV.WeightKey      = 'FLASCALE'; %

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% prepare list of images
IsTmp = false;
if (isstruct(ImList) || SIM.issim(ImList))
    % save SIM as a temp FITS files
    ImList = sim2fits(ImList,'TmpName',true);
    IsTmp  = true;
end
[ListFile,ListCell] = Util.files.create_list(ImList,[]);
Nim = length(ListCell);


% construct ...,key,val, parameters name
ProgFullPath = Util.string.construct_fullpath(InPar.SWarpProg,InPar.SWarpDir,[],mfilename);

if (~isempty(InPar.RA) && ~isempty(InPar.Dec))
    RA  = celestial.coo.convertdms(InPar.RA,'gH','SH');
    Dec = celestial.coo.convertdms(InPar.Dec,'gD','SD');
    CooString = sprintf('-CENTER %s,%s',RA,Dec);
else
    CooString = '';
end
    

KVstr = construct_keyval_string(InPar.SWarpPars);
if (~isempty(InPar.Weight))
    % update FITS header keywords for FLUX scale
    for Iim=1:1:Nim
        fits_write_keywords(ListCell{Iim},{InPar.WeightKey,InPar.Weight(Iim)});
        KVstr = sprintf('%s -FSCALE_KEYWORD %s',KVstr,InPar.WeightKey);  % add flux keyword
    end    
end


if (isempty(InPar.ConfigFile)==1)
   KVstr = sprintf('%s %s',KVstr,CooString);
else
   KVstr = sprintf('-c %s %s %s',InPar.ConfigFile,KVstr,CooString);
end

ExecStr = sprintf('%s @%s %s',ProgFullPath,ListFile,KVstr);
system(ExecStr);


% clean tmp files
delete(ListFile);
if (IsTmp)
    Util.files.delete_cell(ListCell);
end

% output file name
if (~isempty(InPar.OutIm))
    movefile('coadd.fits',InPar.OutIm);
end
if (~isempty(InPar.OutWeightIm))
    movefile('coadd.weight.fits',InPar.OutWeightIm);
end


