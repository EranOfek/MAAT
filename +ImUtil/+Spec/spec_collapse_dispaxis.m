function Collapse=spec_collapse_dispaxis(Sim,varargin)
%--------------------------------------------------------------------------
% spec_collapse_dispaxis function                                   ImSpec
% Description: Given an image/s collapse each image to 1 D vector
%              along one of the dimensions (e.g., "dispersion direction").
% Input  : - Input images to collapse.
%            The following inputs are possible:
%            (1) Cell array of image names in string format.
%            (2) String containing wild cards (see create_list.m for
%                option). E.g., 'lred00[15-28].fits' or 'lred001*.fits'.
%            (3) Structure array of images (SIM).
%                The image should be stored in the 'Im' field.
%                This may contains also mask image (in the 'Mask' field),
%                and an error image (in the 'ErrIm' field).
%            (4) Cell array of matrices.
%            (5) A file contains a list of image (e.g., '@list').
%          * Arbitrary number of ...,key,val,... input arguments.
%            The following keywords are available:
%            'DispDir' - Direction in which to collapse the image {'x'|'y'}
%                        (e.g., dispersion direction). Default is 'x'.
%            'Range'   - Range of image along the dispersion axis [min max]
%                        to collapse. Default is empty. If empty then use
%                        the entire image.
%            'CollapseAlgo' - Collapse algorithm:
%                        'optimal' - Mean weighted by inverse variance
%                                    (including readnoise). Default.
%                        'mean'    - Mean.
%                        'median'  - Median.
%                        'quantile'- Quantile with probability given by
%                                    Prct.
%            'Prct'    - Fraction for the quantile calculation.
%                        Default is 0.9. This may be good for selecting
%                        sources with strong emission lines.
%            'ImageField' - The field name in the structure array
%                        containing the image. Default is '.Im'.
%            'FitsReadOpt' - Cell array of additional parameters to pass
%                        to the fitsread.m function. Default is {}.
%            --- Additional parameters
%            Any additional key,val, that are recognized by images2sim.m.
% Output : - Cell array of collapsed vectors.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Collapse=spec_collapse_dispaxis(rand(1000,500));
%          Collapse=spec_collapse_dispaxis('*.fits','RN',6,'DispDir','y');
%          Collapse=spec_collapse_dispaxis('red0021.fits')
% Reliable: 2
%--------------------------------------------------------------------------


ImageField  = 'Im';
HeaderField = 'Header';
%FileField   = 'ImageFileName';
MaskField   = 'Mask';
BackImField = 'BackIm';
ErrImField  = 'ErrIm';


DefV.ImageField   = ImageField;
DefV.DispDir      = 'x';
DefV.FitsReadOpt  = {};
DefV.Range        = [];  % range to collapse
DefV.CollapseAlgo = 'optimal'; % {'optimal'|'mean'|'median','quantile'}
DefV.RN           = 7;   % [e-]
DefV.Prct         = 0.9;
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

ImageField = InPar.ImageField;

% set dispersion direction
DispDim.x = 2;
DispDim.y = 1;
DispDir = DispDim.(InPar.DispDir);

% set up input files/images
Sim = images2sim(Sim);
Nim = numel(Sim);

Collapse = cell(Nim,1);
for Iim=1:1:Nim,
    % rotate
    switch lower(InPar.DispDir)
        case 'x'
            CurrImage = Sim(Iim).(ImageField);
        case 'y'
            CurrImage = Sim(Iim).(ImageField).';
        otherwise
            error('Unknown DispDir option');
    end
    
    % select sub image
    if (~isempty(InPar.Range)),
        CurrImage = CurrImage(:,InPar.Range(1):InPar.Range(2));
    else
        % do nothing
    end
    
    % collapse image along the dispersion axis
    switch lower(InPar.CollapseAlgo)
        case 'median'
            Collapse{Iim} = median(CurrImage,2);
        case 'mean'
            Collapse{Iim} = mean(CurrImage,2);
        case 'prctile'
            Collapse{Iim} = quantile(CurrImage,InPar.Prct,2);
        case 'optimal'
            Var = CurrImage + InPar.RN.^2;
            Collapse{Iim} = nansum(CurrImage./Var,2)./nansum(1./Var,2);
        otherwise
            error('Unknown CollapseAlgo option');
    end
    
end


    

