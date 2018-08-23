function [R]=coadd_proper(Sim,varargin)
% Proper coaddition of images in a SIM object
% Package: @SIM
% Description: Proper coaddition of images in a SIM object.
%              The following steps will be applied: gain correction,
%              source extraction, registration, estimate variance, estimate
%              zero points, estimate PSFs, subtract background, and apply
%              proper coaddition (Zackay & Ofek 2015b).
% Input  : - A multi element SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Align' - A logical flag indicating if to attempt to align the
%                      images prior to coaddition. Default is true.
%            'AlignPar' - A cell array of additional arguments to pass to
%                      SIM/align.m. Default is {}.
%            'ExecField' - Field to coadd. Default is 'Im'.
%            'GainCorrect' - A logical flag indicating if to correct the
%                      gain to 1. Default is true.
%            'BackReCalc' - A logical flag indicating if to re-calculate
%                      the background and error images. Default is false.
%            'BackPar' - A cell array of additional arguments to pass to
%                      SIM/background.m. Default is {}.
%            'FunGlobalBack' - A function handle to use in order to convert
%                      the local background to global background.
%                      Default is @median.
%            'WeightVar' - A vector of variance weights for each image in
%                      the SIM array to coadd. If empty, then will estimate
%                      the variance from the Err^2 image.
%                      Default is empty.
%            'WeightTran' - A vector of trasperncy weight for each image in
%                      the SIM array to coadd. If empty then will estimate
%                      it using SIM/phot_zp.m.
%                      Default is empty.
%            'PhotZPPar' - A cell array of additional arguments to pass to
%                      SIM/phot_zp.m. Default is {}.
%            'FilterPar' - A cell array of additional arguments to pass to
%                      SIM/filter.m. Default is {}.
%            'CoaddPar' - A cell array of additional arguments to pass to
%                      SIM/coadd.m. Default is {}.
% Output : - A SIM object containing the proper coaddition of the images.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: % Coadd with all defaults
%          R = coadd_proper(S);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField   = SIM.ImageField;   % 'Im';
BackField    = SIM.BackField;
ErrField     = SIM.ErrField;     % 'ErrIm';
WeightField  = SIM.WeightField;  % 'WeightIm';
MaskField    = SIM.MaskField;    % 'Mask';
%MaskDicField = 'MaskDic';

DefV.Align              = true;
DefV.AlignPar           = {};
DefV.ExecField          = SIM.ImageField;
DefV.GainCorrect        = true;
DefV.BackReCalc         = false;
DefV.BackPar            = {'Block',[256 256]};
DefV.FunGlobalBack      = @median;   % Fun defined on SIM and of the form: median(Sim,'ExecField',ErrField)
DefV.WeightVar          = [];
DefV.WeightTran         = [];
DefV.PhotZPPar          = {};
DefV.FilterPar          = {};
DefV.CoaddPar           = {};

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Nsim  = numel(Sim);

R = SIM;
    

%------------------------
%--- correct for gain ---
%------------------------
Sim = gain_correct(Sim);

%------------------------
%--- align the images ---
%------------------------
if (InPar.Align)
    Sim = align(Sim,[],InPar.AlignPar{:});
end
   
% assume the images are registered

%------------------------------
%--- Estimate image weights ---
%------------------------------
% check if background and std exist - if not then calculate
if (InPar.BackReCalc)
    BackExist = false(Nsim,1);
else
    BackExist = isfield_populated(Sim,BackField);
end
Sim(~BackExist) = background(Sim(~BackExist),InPar.BackPar{:});

% global variance of image
GlobalBack = InPar.FunGlobalBack(Sim,'ExecField',BackField);
if (isempty(InPar.WeightVar))
    GlobalVar  = InPar.FunGlobalBack(Sim,'ExecField',ErrField).^2;
else
    % user supplied vector of variances (per image)
    GlobalVar = InPar.WeightVar;
end

% Estimate images transperencies
if isempty(InPar.WeightTran)
    [Res,TranZP] = phot_zp(Sim,InPar.PhotZPPar{:});
    TranZP;
else
    TranZP = InPar.WeightTran;
end

% Sum the original background of the images
% FFU

%---------------------------
%--- Subtract background ---
%---------------------------
SimB = sub_background(Sim);

% Filter image with its PSF
SimB = filter(SimB,InPar.FilterPar{:});

% FFT the images
FFT_SimB = fft2(SimB);

% Shift and pad the PSF such it is in the corner of the full image size
% SP is a SIM array in which the PSF is in the 'Im' field.
ImageSize = fliplr(imagesize(Sim));  % [Y X] imagesize
SP = psf_pad2sim(Sim,ImageSize,true);



% coadd the images using proper coaddition
FFT_M = fft2(SimB);   % fft2 of the measured image
FFT_P = fft2(SP);     % conj(fft2) of the PSF

% Note that the SIM.*vector operation multiply each SIM element by
% the corresponding vector element
FFT_S  = coadd(FFT_M.*conj(FFT_P).*(TranZP(:)./GlobalVar(:)),'MethodZero','none','MethodScale','none');
FFT_PR = sqrt(coadd(FFT_P.*conj(FFT_P).*(TranZP(:).^2./GlobalVar(:)),InPar.CoaddPar{:},'MethodZero','none','MethodScale','none'));
FFT_R  = FFT_S./FFT_PR;

R      = (ifft2(FFT_R));
% populate the PSF P_R in R:
% Note that the PSF is cuurently in the Im field of FFT_PR and
% it has the full image size
R      = setpsf(R,fftshift(ifft2(FFT_PR.(ImageField))));

% FFU
R.BackIm   = [];
R.ErrIm    = [];
R.WeightIm = [];

% [Mom,M2]=ImUtil.Im.im_moments(R.PSF,1025.5 ,2049,3,1.5,3)
%ImageSizeR = size(R.(ImageField));
%CenterPSF = ImageSizeR.*0.5 + 1;
% FFU - Need a could that optionaly shrink the size of the PSF...


%%% NEED TO DEAL with back, err and MASK

