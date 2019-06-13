function [Mask,CleanImage]=crdetect_lacos(Image,varargin)
% L.A. cosmic cosmic ray detection for astronomical images
% Package: ImUtil.Im
% Description: Find and remove cosmic rays in an astronomical image
%              using the L.A.cosmic algorithm (van Dokkum 2001).
% Input  : - Image matrix.
%          * Arbitrary number of pairs of input arguments ...,key,val,...
%            The following keywords are available:
%            'RepInf'- Replace infinities with a constant.
%                     If empty, do nothing, otherwise this is the constant.
%                     Default is 1e6.
%            'Gain' - CCD gain for noise estimation. Default is 1.
%            'RN'   - CCD read noise [e-] for noise estimation.
%                     Default is 10.
%                     This can be an image of readnoise in e- per pixel.
%            'Nsigma' - CR detection threshold in sigma. Default is 10.
%            'Fth'    - Fine structure threshold. Default is 2.
%                       If FWHM is provided than this parameter is
%                       overrided.
%            'FWHM'   - PSF FWHM to estimate 'Fth' based on Figure 4
%                       in van Dokkum (2001).
%            'Mask'   - Bit mask in which to set a bit for each pixel
%                       which is affected by CR.
%                       Default is empty. If empty then will create
%                       a new mask.
%            'Bit_CR' - Index of bit in the bit mask to flag.
%                       Either an index or a function handle
%                       that get 'Bit_CR' and returns the Bit_CR index.
%                       If empty then will return a logical mask
%                       (true for CR, false for no CR).
%                       Default is empty.
%            'BWmorph'- If not empty, then increase the CR effected area
%                       by running a morphological filter using bwmorph.m.
%                       Default is 'dilate'.
%                       Other useful options are: 'majority' and 'bridge'.
%            'BWmorphN'-Number of times to run the morphological filter.
%            'MaskType' - Type of bit mask to create. Default is 'uint16'.
%            'IntMethod'- inpaint_nans.m interpolation method. Default is 2.
% Output : - Bit mask of flag image.
%            If 'Bit_CR' is not provided than this is a flag image (logicals)
%            in which pixels affected by CR are true, all the rest are
%            false.
%            If 'Bit' is not empty, then this is a bit mask in which
%            of the CR affected pixels.
%          - Image interpolated over the CRs.
% Reference: http://www.astro.yale.edu/dokkum/lacosmic/
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: A=rand(6,6); A(3,3)=1000;
%          [Mask,CleanImage]=ImUtil.Im.crdetect_lacos(A);
% Reliable: 2
%--------------------------------------------------------------------------
import Util.array.*
import Util.stat.*

DefV.RepInf   = 1e6;
DefV.Gain     = 1;
DefV.RN       = 10;
DefV.Nsigma   = 8;
DefV.Fth      = 2;
DefV.FWHM     = []; %2.5;  % pix
DefV.Mask     = [];
DefV.Bit_CR   = [];
DefV.BWmorph  = 'dilate';
DefV.BWmorphN = 1;
DefV.MaskType = 'uint16';
DefV.IntMethod= 2;
DefV.NoiseMethod = 'poisson'; % {'poisson','mode_fit'}
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Ftable = [0.9 10;1 6;1.5 3.5; 2 2; 2.5 1.5; 3 1; 5 0.9; 10 0.8];
if (~isempty(InPar.FWHM))
    InPar.Fth = interp1(Ftable(:,1),Ftable(:,2),InPar.FWHM,'linear');
end

switch lower(InPar.NoiseMethod)
    case 'poisson'
        NoiseImage = sqrt(InPar.Gain.*abs(Image) + InPar.RN.^2)./sqrt(InPar.Gain);
    case 'mode_fit'
        [~,NoiseImage] = mode_fit(Image);
    otherwise
        error('Unknown NoiseMethod option');
end

% replace inf
if (~isempty(InPar.RepInf))
    Image(isinf(Image)) = InPar.RepInf;
end  

% subsample
Image2 = imresize(Image,2,'nearest');
%Image2 = imresize(Image,2,'bicubic');

% imlaplacian
Lap2 = ImUtil.Im.imlaplacian(Image2);

% replace L2<0 with 0
Lap2(Lap2<0) = 0;

% convolve with [1 1;1 1]

% return to original sampling
%LAC = imresize(Lap2,0.5,'nearest');
LAC = imresize(Lap2,0.5,'bicubic');
%LAC = Lap2(2:2:end,2:2:end);

S = LAC./(NoiseImage.*sqrt(5./4)); % .*Nsigma;
St = S - medfilt2(S,[5 5]);

% fine structure image
F3 = medfilt2(Image,[3 3]);
F  = F3 - medfilt2(F3,[7 7]);
ImageCR = LAC./F>InPar.Fth & St>InPar.Nsigma;

if (~isempty(InPar.BWmorph))
    ImageCR = bwmorph(ImageCR,InPar.BWmorph,InPar.BWmorphN);
    %ImageCR = imdilate(ImageCR,ones(InPar.Increase));
    %ImageCR(medfilt2(ImageCR,[InPar.Increase])>0) = true;
end    

% bit index
InPar.Bit_CR = get_bitmask_def(InPar.Bit_CR,'Bit_CR_LACOS');

if (isempty(InPar.Bit_CR))
    % Mask output is a logical image
    Mask = ImageCR;
else
    % Mask image is a bit mask
    if (isempty(InPar.Mask))
        % create Mask
        Mask = maskflag_set([],InPar.MaskType,InPar.Bit_CR,ImageCR);
    else
        % set Mask
        Mask = maskflag_set(InPar.Mask,InPar.MaskType,InPar.Bit_CR,ImageCR);
    end
end

if (nargout>1)
   CleanImage = ImUtil.Im.iminterp(Image,ImageCR,'IntMethod',InPar.IntMethod);
end


