function spec_extract_1d(Sim,varargin)
%--------------------------------------------------------------------------
% spec_extract_1d function                                          ImSpec
% Description: Given a background subtracted 2D image in which the
%              spectrum is aligned along one axis, extract (or fit)
%              the 1D spectrum.
% Input  : - A single image in one of the following forms:
%            (1) A structure array (SIM).
%            (2) A file name.
%            (3) A matrix.
%            (4) A cell array with a single file name string.
%            (5) A cell array with a single matrix image.
%            See image2sim.m for options.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:

%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
%            image2sim.m


import Util.array.*

ImageField  = 'Im';
HeaderField = 'Header';
FileField   = 'ImageFileName';
MaskField   = 'Mask';
%BackImField = 'BackIm';
%ErrImField  = 'ErrIm';

DefV.DispDir    = 'x';
DefV.SpecCenter = [];
DefV.MethodExt1 = 'aper';   % {'center'|'aper'|'psf'|'gauss'}
DefV.AperRad    = 3;   % aperture radius or PSF radius
DefV.PSF        = [];
DefV.RangePSF   = [];


DefV.Gain       = 1;
DefV.RN         = 10;  % [e-]

InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

% read image
Sim = image2sim(Sim,varargin{:});

% make sure dispersion axis is along the x-axis
switch lower(InPar.DispDir)
    case 'y'
        Sim.(ImageField) = Sim.(ImageField).';
    otherwise
        % do nothing
end

% image size
SizeIm = size(Sim.(ImageField));
Ndisp  = SizeIm(2);
if (isempty(InPar.RangePSF))
    InPar.RangePSF = [1 SizeIm(2)];
end

% spectrum center
if (isempty(InPar.SpacCenter)),
    % assume the spectrum is in the image center
    InPar.SpecCenter = (SizeIm(1) + 1).*0.5;
end

% pixel indices
ExtractSpec.Pix = (1:1:SizeIm(2));

% extract spectrum
switch lower(InPar.MethodExt1)
    case 'center'
        % extract central line (in DN)
        ExtractSpec.Spec = Sim.(ImageField)(InPar.SpacCenter,:);
        
        if (isfield(Sim,MaskField)),
            ExtractSpec.Mask = Sim.(MaskField)(InPar.SpacCenter,:);
        end
        if (isfield(Sim,BackImField)),
            ExtractSpec.Back = Sim.(BackImField)(InPar.SpacCenter,:);
        else
            ExtractSpec.Back = zeros(size(ExtractSpec.Spec));
        end
        
        % estimate error (in DN)
        ExtractSpec.SpecErr =  sqrt(ExtractSpec.Spec.*InPar.Gain + ExtractSpecBack.*InPar.Gain + InPar.RN.^2)./InPar.Gain;
        
    case 'aper'
        AperInd = (floor(InPar.SpacCenter-InPar.AperRad):ceil(InPar.SpacCenter+InPar.AperRad));
        ExtractSpec.Spec = nansum(Sim.(ImageField)(AperInd,:));
        
        if (isfield(Sim,MaskField)),
            ExtractSpec.Mask = sum_bitor(Sim.(MaskField)(AperInd,:),1);
        end

        if (isfield(Sim,BackImField)),
            ExtractSpec.Back = sum(Sim.(BackImField)(AperInd,:));
        else
            ExtractSpec.Back = zeros(size(ExtractSpec.Spec));
        end
        
        % estimate error (in DN)
        ExtractSpec.SpecErr =  sqrt(ExtractSpec.Spec.*InPar.Gain + ExtractSpecBack.*InPar.Gain + InPar.RN.^2)./InPar.Gain;
        
    case 'psf'
        AperInd = (floor(InPar.SpacCenter-InPar.AperRad):ceil(InPar.SpacCenter+InPar.AperRad));

        % construct PSF
        if (isempty(InPar.PSF)),
            % construct PSF by collapsing the spectrum
            PSF = nanmedian(Sim.(ImageField)(AperInd,InPar.RangePSF(1):InPar.RangePSF(2)),2);
        else
            PSF = InPar.PSF;
        end
        % normalize PSF
        PSF = PSF./sum(PSF);
        
        % fit the PSF
        ExtractSpec.Spec = zeros(1,Ndisp);
        for Idisp=1:1:Ndisp,
            ExtractSpec.Spec(Idisp) = Sim.(ImageField)(AperInd,Idisp)./PSF;
        end
        
    otherwise
        error('Unknown MethodExt1 option');
end
    