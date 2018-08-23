function [Sim,Npixsat]=flag_saturated(Sim,varargin)
% Flag saturated pixels in bit mask image
% Package: @SIM
% Description: Get saturation level for image and set bit mask of
%              saturated pixels.
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ImageField' - The image field in the SIM object in which to
%                           look for saturated pixels. Default is 'Im'.
%            'SatKey'     - A cell array of header keywords in which to
%                           search for the saturation level of the image.
%                           See getkey_fromlist.m for details.
%                           Default is
%                           {'SATURVAL','SATLEVEL','SATURAT','SATUR','SATURATE'}.
%            'SatLevel'   - Default saturation level to use if the
%                           saturation level header keyword doesn't exist.
%                           Default is 50000.
%            'MaskDic'    - Mask dictionary.
%                           Default is @MASK.def_bitmask_pipeline.
%            'Bit_ImSaturated'- Saturation bit mask name in the dictionary
%                           or its index. Default is 'Bit_ImSaturated'.
%            'BitClass'   - Default bit mask class. Default is 'uint32'.
%            'SetNaN      - Set the value of saturated pixels in the 'Im'
%                           field to NaN.
% Output : - A SIM object with the MASK field updated regarding saturated
%            pixels.
%          - An array in which each element indicate the number of
%            saturated pixels in the image.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [S1,Npixsat]=flag_saturated(S);
% Reliable: 2
%--------------------------------------------------------------------------


ImageField       = 'Im';
MaskDicField     = 'MaskDic';

DefV.ImageField         = ImageField;
DefV.SatKey             = {'SATURVAL','SATLEVEL','SATURAT','SATUR','SATURATE'};
DefV.SatLevel           = 50000;
DefV.MaskDic            = @MASK.def_bitmask_pipeline;
DefV.Bit_ImSaturated    = 'Bit_ImSaturated';
DefV.BitClass           = 'uint32';
DefV.SatNaN             = false;

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);



% get stauration level:
SatLevel = cell2mat(getkey_fromlist(Sim,InPar.SatKey));

% Deal with NaN (Saturation level is not in header)
if (any(isnan(SatLevel)))
    warning('Some saturation level are NaN - use default value of %f instead',InPar.SatLevel);
    SatLevel(isnan(SatLevel)) = InPar.SatLevel;
end

Nsim = numel(Sim);

Npixsat = nan(size(Sim));
for Isim=1:1:Nsim
    % Set Mask dictionary
    if (isempty(Sim(Isim).(MaskDicField)))
        Sim(Isim).(MaskDicField) = InPar.MaskDic;
    end
    
    % Search for saturated pixels
    Flag = Sim(Isim).(ImageField) > SatLevel(Isim);
    
    % set bit maks
    Sim(Isim) = bitmask_set(Sim(Isim),Flag,InPar.Bit_ImSaturated,InPar.BitClass,@bitor);
    
    % count saturated pixels
    if (nargout>1)
        Npixsat(Isim) = sum(Flag(:));
    end
    
    % set pixel value to NaN
    if (InPar.SatNaN)
        Sim(Isim).(ImageField)(Flag) = NaN;
    end
end

    
    
    
    