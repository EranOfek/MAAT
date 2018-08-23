function [Out,Norm]=spec_coadd1d(InSpec,varargin)
%--------------------------------------------------------------------------
% spec_coadd1d function                                             ImSpec
% Description: 
% Input  : - A structure array that contains a spectra in each cell.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Coadd'     - One of the following coaddition methods:
%                          'median' - median (default).
%                          'mean'   - mean.
%                          'wmean'  - mean weighted by the errors.
%                          Note that the mask will always be coadd using
%                          or operation.
%            'NormCoadd' - Normalize spectra before coaddition.
%                          {'none','lsq'}. Default is 'lsq'.
%            'FluxField' - The field name in the input structure array
%                          that contains the flux. Default is 'Flux'.
%            'WaveField' - The field name in the input structure array
%                          that contains the wavelength. Default is 'Wave'.
%            'ErrField'  - The field name in the input structure array
%                          that contains the wavelength. Default is 'FluxErr'.
%            'MaskField' - The field name in the input structure array
%                          that contains the mask. Default is 'Mask'.
%            'Wave'      - Vector of wavelengths to which to interpolate
%                          all the spectra. If empty, then use first
%                          input spectra. Default is empty.
%            'InterpMethod' - Interpolation method. Default is 'linear'.
%                          This interpolation method is used for
%                          interpolating the flux and error, but the mask
%                          is interpolated using 'nearest'.
%            'MaskType'  - Mask type. Default is 'uint16'.
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - A structure containing the following fields:
%            .Wave - wavelength
%            .Flux
%            .Mask
%            .ErrStd
%            .ErrStdN
%            .ErrW
%          - Vector of normalization that was used for each spectra.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Out=spec_coadd1d(InSpec,'WaveField','PrimWC');
% Reliable: 2
%--------------------------------------------------------------------------
import Util.array.*

DefV.FluxField    = 'Flux';
DefV.ErrField     = 'FluxErr';
DefV.MaskField    = 'Mask';
DefV.WaveField    = 'Wave';
DefV.Wave         = [];
DefV.InterpMethod = 'linear';
DefV.MaskType     = 'uint16';
DefV.Coadd        = 'median';
DefV.NormCoadd    = 'lsq';

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

%if (iscell(InSpec)),
%    InSpec = cell2struct(InSpec,InPar.SpecField,1);
%end
%Col = cell2struct(num2cell(1:1:length(InPar.Columns)),InPar.Columns,2);

if (isempty(InPar.Wave))
    InPar.Wave = InSpec(1).(InPar.WaveField);
end
Nwave = numel(InPar.Wave);
Wave  = InPar.Wave;

Nspec     = numel(InSpec);
Mat.Flux  = zeros(Nwave,Nspec).*NaN;
Mat.Err   = zeros(Nwave,Nspec).*NaN;
Mat.Mask  = zeros(Nwave,Nspec,InPar.MaskType);

for Ispec=1:1:Nspec
    Mat.Flux(:,Ispec) = interp1(InSpec(Ispec).(InPar.WaveField), InSpec(Ispec).(InPar.FluxField), Wave,InPar.InterpMethod);

    if (isfield(InSpec,InPar.ErrField))
        Mat.Err(:,Ispec)  = interp1(InSpec(Ispec).(InPar.WaveField), InSpec(Ispec).(InPar.ErrField), Wave,InPar.InterpMethod);
    end
    if (isfield(InSpec,InPar.MaskField))
        Mat.Mask(:,Ispec) = interp1(InSpec(Ispec).(InPar.WaveField), double(InSpec(Ispec).(InPar.MaskField)), Wave,'nearest');
        Mat.Mask(:,Ispec) = cast(Mat.Mask(:,Ispec),InPar.MaskType);
    end
end

% coadd data
if (all(isnan(Mat.Err)))
    Err = ones(size(Mat.Err));
else
    Err = Mat.Err;
end

% Normalize prior to coaddition
switch lower(InPar.NormCoadd)
    case {'no','none'}
        % do nothing
    case 'lsq'
        if (Nspec==1)
            Norm = 1;
        else
            IndNN = ~isnan(sum(Mat.Flux,2));
            H     = Mat.Flux(IndNN,:);
            Norm  = H(:,2:end)\H(:,1);
            Norm  = [1; Norm].';
        end
        Mat.Flux = bsxfun(@times,Mat.Flux,Norm);
        Mat.Err  = bsxfun(@times,Mat.Err,Norm);
       
    otherwise
        error('Unknown NormCoadd method');
end


switch lower(InPar.Coadd)
    case 'median'
        Flux = nanmedian(Mat.Flux,2);
    case 'mean'
        Flux = nanmean(Mat.Flux,2);
    case 'wmean'
        % weighted mean
        Flux = nansum(Mat.Flux./(Err.^2),2)./nansum(1./(Err.^2),2);
    otherwise
        error('Unknown Coadd option');
end
Out.Wave    = Wave;
Out.Flux    = Flux;
Out.ErrStd  = nanstd(Err,[],2);
Out.ErrStdN = nanstd(Err,[],2)./sqrt(sum(~isnan(Mat.Flux),2));
Out.ErrW    = sqrt(nansum(1./(Err.^2),2));

Out.Mask    = sum_bitor(Mat.Mask,2);
        
        

    
    
    