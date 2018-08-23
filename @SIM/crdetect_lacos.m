function Sim=crdetect_lacos(Sim,varargin)
% Detect cosmic rays using the Laplasian Edge Detector (van Dokkum 2001)
% Package: SIM
% Description: Find and remove cosmic rays in an astronomical image
%              using the L.A.cosmic algorithm (van Dokkum 2001).
% Input  : - SIM object.
%          * Arbitrary number of pairs of input arguments ...,key,val,...
%            The following keywords are available:
%            'InterpOverCR' - Interpolate the image over CRs.
%                           Default is true.
%            'CR_BitName'   - Bit mask name to switch on in case of a CR.
%                           Default is 'Bit_CR_LACOS'.
%            'BackPar'      - Cell array of additional arguments to pass to
%                           SIM/background.
%            'GainCorrect'  - Correct image to gain=1. Default is true.
%            'GainCorrectPar' - Cell array of additional arguments to pass
%                           to SIM/gain_correct.
%            'ImageField'   - SIM object image field on which to operate.
%                           Default is 'Im'.
%            'RepInf'- Replace infinities with a constant.
%                     If empty, do nothing, otherwise this is the constant.
%                     Default is 1e6.
%            'Nsigma' - CR detection threshold in sigma. Default is 10.
%            'Fth'    - Fine structure threshold. Default is 2.
%                       If FWHM is provided than this parameter is
%                       overrided.
%            'FWHM'   - PSF FWHM to estimate 'Fth' based on Figure 4
%                       in van Dokkum (2001).
%            'BWmorph'- If not empty, then increase the CR effected area
%                       by running a morphological filter using bwmorph.m.
%                       Default is 'dilate'.
%                       Other useful options are: 'majority' and 'bridge'.
%            'BWmorphN'-Number of times to run the morphological filter.
%            'IntMethod'- inpaint_nans.m interpolation method. Default is 2.
% Output : - A SIM object in which the MASK image is updated with candidate
%            cosmic ray locations, and optionaly, the image is interpolated
%            over these locations.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=FITS.read2sim('*.fits');
%          S=crdetect_lacos(S);
% Reliable: 2

import Util.array.*
import Util.stat.*

ErrField      = SIM.ErrField;
ImageField    = SIM.ImageField;

DefV.InterpOverCR       = true;
DefV.CR_BitName         = 'Bit_CR_LACOS';
DefV.BackPar            = {};
DefV.GainCorrect        = true;
DefV.GainCorrectPar     = {};
DefV.ImageField         = SIM.ImageField;
DefV.RepInf             = 1e6;
DefV.Nsigma             = 10;
DefV.Fth                = 2;
DefV.FWHM               = 6; %[]; %2.5;  % pix
DefV.Mask               = [];
DefV.BWmorph            = 'dilate';
DefV.BWmorphN           = 1;
DefV.IntMethod          = 2;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% check if noise image exist
IsErrExist = isfield_populated(Sim,ErrField);
Sim(~IsErrExist) = background(Sim(~IsErrExist),InPar.BackPar{:});


Ftable = [0.9 10;1 6;1.5 3.5; 2 2; 2.5 1.5; 3 1; 5 0.9; 10 0.8];
if (~isempty(InPar.FWHM))
    InPar.Fth = interp1(Ftable(:,1),Ftable(:,2),InPar.FWHM,'linear');
end

% correct Gain
Sim = gain_correct(Sim,InPar.GainCorrectPar{:});



Nsim = numel(Sim);

for Isim=1:1:Nsim
    % for each image
    Image = Sim(Isim).(InPar.ImageField);
    
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

    S = LAC./(Sim(Isim).(ErrField).*sqrt(5./4)); % .*Nsigma;
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

    % update mask image
    Sim(Isim) = bitmask_set(Sim(Isim),ImageCR,InPar.CR_BitName); %,BitType,CombFun)
    
    % interpolate over CRs
    if (InPar.InterpOverCR)
        Sim(Isim).(ImageField) = ImUtil.Im.iminterp(Image,ImageCR,'IntMethod',InPar.IntMethod);
    end
    
end



