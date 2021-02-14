function [SN]=sn_spec(Spec,varargin)
% S/N calculator for long-slit spectra
% Package: telescope.sn
% Description: Simulate long-slit spectral observations and estimate the
%              S/N per resolution element.
% Input  : - Spectrum [Wave(Ang), Flux(cgs/A)] or an AstSpec object
%            containing spectrum. Default is 'QSO_SDSS' at mag 20.0
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'RN'      - Readout noise [e-]. Default is 10.
%            'Gain'    - CCD gain [e-/ADU] (for digitization noise).
%                        Default is 1.5.
%            'ExpTime' - Total exposure time [s]. Default is 600.
%            'SubExp'  - Number of sub exposure within total exposure time.
%                        Default is 1.
%            'DC'      - CCD dark curent [e-/pix/s]. Default is 1e-3.
%            'Back'    - A background spectrum [Wave, Flux(cgs/A/arcsec^2)]
%                        Either an AstSpec object a matrix or background
%                        catalog name in the cats.spec.SkyBack pacakge.
%                        Default is 'Gemini_SkyBack_dark'.
%            'PixScale' -Pixel scale in spatial direction [arcsec/pix].
%                        Default is 0.6.
%            'Seeing'   -Seeing FWHM [arcsec]. Default is 1.2.
%            'SpecRes'  - 
%            'SpecSampling' -
%            'AirMass'  - 
%            'AtmosphericExt'
%            'Aper'
%            'TelTh'
%            'SpecTh'
%            'QE'
%            'ExtractionSemiWidth'
%            'SpecRange'
%            'SpecFluxUnits'        = 'cgs/A';
%            'BackFluxUnits'        = 'cgs/A';
%            'SpecWaveUnits'        = 'A';
%            'BackWaveUnits'        = 'A';
%            'ColW'                 = 1;
%            'ColF'                 = 2;
%            'InterpMethod'         = 'linear';
% Output : - A structure array with the following fields:
%            'Wave' - A vector of wavelength (Ang) at which the S/N is
%                     calculated. This vector is is given by the spectral
%                     sampling.
%            'SNperResEl' - S/N per resolution element as extracted using
%                     optimal extraction. Note that this is not the Horne
%                     et al. extraction.
%            'ResEl' - Resolution element (Ang).
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Spec = AstSpec.get_galspec('Gal_E');
%          Spec = scale2mag(Spec,20.5);
%          Spec.synphot('SDSS','r','AB')
%          SN=telescope.sn.sn_spec(Spec)
% Reliable: Under tests
%--------------------------------------------------------------------------

if (nargin==0)
    Spec = [];
end
if (isempty(Spec))
    Spec = AstSpec.get_galspec('QSO_SDSS');
    Spec = scale2mag(Spec,20);
end

DefV.RN                   = 3;   % [e-]
DefV.Gain                 = 1.5;  % [e-/ADU]
DefV.ExpTime              = 900;  % [s]
DefV.SubExp               = 1;
DefV.SlitWidth            = 3;    % [arcsec]
DefV.ExtraEff             = 0.8;  % Additional inefficiencies efficiency
DefV.DC                   = 1e-2; % [e-/pix/s]
DefV.Back                 = 'Gemini_SkyBack_dark';   % [cgs/A per arcsec]
DefV.PixScale             = 0.6;  % [arcsec/pix]
DefV.Seeing               = 2;    % FWHM [arcsec]
DefV.SpecRes              = @(X) X./3000; % @(X) 10+0.5.*X./200;  % 25;    % FWHM [Ang], or vec or fun
DefV.SpecSampling         = 0.8;      % Ang/pix; [] don't resample
DefV.AirMass              = 1.4;   % 0 for space
DefV.AtmosphericExt       = 'KPNO';
DefV.Aper                 = 254;    % cm
DefV.TelTh                = 0.85;   % or [Wave, Th]
DefV.SpecTh               = 0.6;   % or [Wave, Th]
DefV.QE                   = [3000 0.2; 3500 0.45; 4000 0.9; 4200 0.92;5100 0.88; 7600 0.93;8500 0.9;9000 0.8;10000 0.3]; %0.9;   % or [Wave, Th]
DefV.ExtractionSemiWidth  = 5;    % spatial pixels in which to fit PSF
DefV.SpecRange            = [3500 9000];
DefV.SpecFluxUnits        = 'cgs/A';
DefV.BackFluxUnits        = 'cgs/A';
DefV.SpecWaveUnits        = 'A';
DefV.BackWaveUnits        = 'A';
DefV.ColW                 = 1;
DefV.ColF                 = 2;
DefV.InterpMethod         = 'linear';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% Slit efficiency
SeeingSigma = InPar.Seeing./2.35;
SlitEff     = 2.*(normcdf(0.5.*InPar.SlitWidth./SeeingSigma,0,1)-normcdf(0,0,1));


% read spectrum
if (AstSpec.isastspec(Spec))
    Spec = astspec2mat(Spec);
end

% resample spectrum at 
if (isempty(InPar.SpecSampling))
    % do not change resampling
else
    StartW = min(Spec(:,InPar.ColW));
    EndW   = max(Spec(:,InPar.ColW));
    
    SamplingVec = (StartW:InPar.SpecSampling:EndW)';
    
    % need to interpolate and sum flux within each resolution element
    SamplingO           = diff(Spec(:,InPar.ColW));
    SamplingO           = [SamplingO(1); SamplingO; SamplingO(end)];
    SamplingO           = 0.5.*(SamplingO(1:end-1) + SamplingO(2:end));

    SamplingN           = diff(SamplingVec);
    SamplingN           = [SamplingN(1); SamplingN; SamplingN(end)];
    SamplingN           = 0.5.*(SamplingN(1:end-1) + SamplingN(2:end));

    SamplingO           = interp1(Spec(:,InPar.ColW), SamplingO, SamplingVec, InPar.InterpMethod);
    
    %Spec        = [SamplingVec, SamplingN./SamplingO .* interp1(Spec(:,InPar.ColW),Spec(:,InPar.ColF),SamplingVec,InPar.InterpMethod)];
    % convert for per pixel later
    Spec        = [SamplingVec, interp1(Spec(:,InPar.ColW),Spec(:,InPar.ColF),SamplingVec,InPar.InterpMethod)];
        
end

% load background spectrum
if (ischar(InPar.Back))
    Back = cats.spec.SkyBack.(InPar.Back);
    Back = astspec2mat(Back); % [erg/cm^2/s/Ang]
elseif (AstSpec.isastspec(InPar.Back))
    Back = astspec2mat(InPar.Back);
else
    % assume a matrix in cgs units
    Back = InPar.Back;  % [erg/cm^2/s/Ang]
end

% resample background spectrum [erg/cm^2/s/Ang]
%Back = [Spec(:,InPar.ColW), SamplingN./SamplingO .* interp1(Back(:,1),Back(:,2),Spec(:,InPar.ColW),InPar.InterpMethod)];
% convert to per pix later...
Back = [Spec(:,InPar.ColW), interp1(Back(:,1),Back(:,2),Spec(:,InPar.ColW),InPar.InterpMethod)];

     
% convert spectrum to ph/A  [ph/cm^2/s/Ang]
Spec(:,InPar.ColF) = convert.flux(Spec(:,InPar.ColF),InPar.SpecFluxUnits,'ph/A',Spec(:,InPar.ColW),InPar.SpecWaveUnits);
Back(:,2)          = convert.flux(Back(:,2),         InPar.BackFluxUnits,'ph/A',Back(:,1),         InPar.BackWaveUnits);
% convert wavelength to ang
% assume spectrum in Ang


% get atmospheric extinction [mag/airmass]
AtmExt = AstSpec.get_atmospheric_extinction(InPar.AtmosphericExt,'mat');
% resample the atmospheric extinction
AtmExt = [Spec(:,InPar.ColW), interp1(AtmExt(:,1),AtmExt(:,2),Spec(:,InPar.ColW),InPar.InterpMethod)];
% apply atmopsheric extinction
Spec(:,InPar.ColF) = Spec(:,InPar.ColF).*  10.^(-0.4.*AtmExt(:,2).*InPar.AirMass);

% telescope aperture [ph/s/Ang]
Spec(:,InPar.ColF) = Spec(:,InPar.ColF).*pi.*(InPar.Aper.*0.5).^2;
Back(:,2)          = Back(:,2)         .*pi.*(InPar.Aper.*0.5).^2;

% telescope throughput [ph/s/Ang]
if (numel(InPar.TelTh)==1)
    % scalar throughput
    Spec(:,InPar.ColF) = Spec(:,InPar.ColF).*InPar.TelTh;
    Back(:,2)          = Back(:,2)         .*InPar.TelTh;
else
    % resample throughput[ph/s/Ang]
    [~,InPar.TelTh] = AstroUtil.spec.eq_sampling(Spec(:,[InPar.ColW, InPar.ColF]),InPar.TelTh,InPar.InterpMethod);
    Spec(:,InPar.ColF) = Spec(:,InPar.ColF).*InPar.TelTh(:,2);
    Back(:,2)          = Back(:,2)         .*InPar.TelTh(:,2);
end

% spectrograph throughput [ph/s/Ang]
if (numel(InPar.SpecTh)==1)
    % scalar throughput
    Spec(:,InPar.ColF) = Spec(:,InPar.ColF).*InPar.SpecTh;
    Back(:,2)          = Back(:,2)         .*InPar.SpecTh;
else
    % resample throughput
    [~,InPar.SpecTh] = AstroUtil.spec.eq_sampling(Spec(:,[InPar.ColW, InPar.ColF]),InPar.SpecTh,InPar.InterpMethod);
    Spec(:,InPar.ColF) = Spec(:,InPar.ColF).*InPar.SpecTh(:,2);
    Back(:,2)          = Back(:,2)         .*InPar.SpecTh(:,2);
end

% CCD QE  [ph/s/Ang]
if (numel(InPar.QE)==1)
    % scalar throughput
    Spec(:,InPar.ColF) = Spec(:,InPar.ColF).*InPar.QE;
    Back(:,2)          = Back(:,2)         .*InPar.QE;
else
    % resample throughput
    %[~,InPar.QE] = AstroUtil.spec.eq_sampling(Spec(:,[InPar.ColW, InPar.ColF]),InPar.QE,InPar.InterpMethod);
    InPar.QE = [Spec(:,InPar.ColW), interp1(InPar.QE(:,1),InPar.QE(:,2),Spec(:,InPar.ColW),InPar.InterpMethod)];
    Spec(:,InPar.ColF) = Spec(:,InPar.ColF).*InPar.QE(:,2);
    Back(:,2)          = Back(:,2)         .*InPar.QE(:,2);
end

% spectrograph sampling
Sampling           = diff(Spec(:,InPar.ColW));
Sampling           = [Sampling(1); Sampling; Sampling(end)];
Sampling           = 0.5.*(Sampling(1:end-1) + Sampling(2:end));   % [Ang/pix]

Spec(:,InPar.ColF) = Spec(:,InPar.ColF).*Sampling;   % [ph/pix/s]
Back(:,2)          = Back(:,2).*InPar.SlitWidth;      % ph in disperion direction per slit width per Ang
Back(:,2)          = Back(:,2)         .*Sampling;   % ph in the disperion dir per slit width and per pix in the disperion dir
SlitWidthPix       = InPar.SlitWidth./InPar.PixScale;   % slit width in pix

% exposure time
% this is the spectrum in e- per spatial-direction pix 
Spec(:,InPar.ColF) = Spec(:,InPar.ColF).*InPar.ExpTime;  % ph/pix
Back(:,2)          = Back(:,2)         .*InPar.ExpTime;  % ph/pix

% convolve the spectrum with the spatial PSF to create a 2-D spectrum
Sigma = InPar.Seeing./2.35;  % FWHM to Gaussian sigma [arcsec]
Sigma = Sigma./InPar.PixScale; % [pixels]
X = (-InPar.ExtractionSemiWidth:1:InPar.ExtractionSemiWidth);
Nx = numel(X);
Y = exp(-X.^2./(2.*Sigma.^2))./(Sigma.*sqrt(2.*pi));
Y = Y./sum(Y);
Spec2D = bsxfun(@times,Spec(:,InPar.ColF),Y);

% convert SpecRes to a vector
if (isa(InPar.SpecRes,'function_handle'))
    % SpecRes is a function handle
    SpecRes = InPar.SpecRes(Spec(:,InPar.ColW));
else
    if (numel(InPar.SpecRes)==1)
        % SpecRes is a scalar
        SpecRes = InPar.SpecRes.*ones(size(Spec(:,InPar.ColW)));
    else
        % assume SpecRes is a vector
        SpecRes = InPar.SpecRes;
    end
end

% convolve with resolution
% Here Spec2D is in units of [ph/pix]
% SpecRes is in units of [Ang]
% Use InPar.SpecSampling to convert Ang to pix
Spec2D = Util.filter.conv1_vargauss(Spec2D,SpecRes./Sampling);

% convolve Back with resolution
Back(:,2) = Util.filter.conv1_vargauss(Back(:,2),SpecRes./Sampling);


% calculate SN per pix

% calculate background variance
BackVar  = Back(:,2) + SlitWidthPix.*InPar.SubExp.*InPar.RN.^2 + SlitWidthPix.*InPar.DC.*InPar.ExpTime + SlitWidthPix.*(InPar.Gain.*0.3).^2;

SNperPix2D = InPar.ExtraEff.*SlitEff*Spec2D./sqrt(bsxfun(@plus,InPar.ExtraEff.*SlitEff.*Spec2D,BackVar));
% SN for detection:
%SNperPix2Dd = InPar.ExtraEff.*SlitEff*Spec2D./sqrt(BackVar);

% add the SN^2...
SN.Wave       = Spec(:,InPar.ColW);
SN.SNperPix   = sqrt(sum(SNperPix2D.^2,2));
SN.ResEl      = SpecRes;
SN.SNperResEl = SN.SNperPix.*sqrt(SN.ResEl./InPar.SpecSampling);




%%
if 1==0
    
    
    VecSpecTh = (0.5:0.025:0.9).';
    VecRN     = 3.*sqrt(1:1:18).';
    Nst = numel(VecSpecTh);
    Nrn = numel(VecRN);
    
    ResMag = nan(Nst,Nrn);
    for Ist=1:1:Nst
        SpecTh = VecSpecTh(Ist);
        for Irn=1:1:Nrn
            RN = VecRN(Irn);
    
            Spec = AstSpec.get_galspec('Gal_E');
            Spec = scale2mag(Spec,20.5);
            SN=telescope.sn.sn_spec(Spec,'RN',RN,'SpecTh',SpecTh);
            F = SN.Wave>5000 & SN.Wave<8000;
            MSN = mean(SN.SNperResEl(F));
            MagDiff = 2.5.*log10((10./MSN));
            NewMag  = 20.5 - MagDiff;
            Spec = scale2mag(Spec,NewMag );
            SN=telescope.sn.sn_spec(Spec,'RN',RN,'SpecTh',SpecTh);
            F = SN.Wave>5000 & SN.Wave<8000;
            MSN = mean(SN.SNperResEl(F));
            MagDiff = 2.5.*log10((10./MSN));
            NewMag = NewMag - MagDiff;
            Spec = scale2mag(Spec,NewMag );
            SN=telescope.sn.sn_spec(Spec,'RN',RN,'SpecTh',SpecTh);
            F = SN.Wave>5000 & SN.Wave<8000;
            MSN = mean(SN.SNperResEl(F));
    
            ResMag(Ist,Irn) = NewMag;
        end
    end
    contour(VecSpecTh,VecRN,ResMag'); shading interp; colorbar 
    set(gca,'XS','log','YS','log')

    H = xlabel('Spectrograph throughput');
    H.Interpreter = 'latex';
    H.FontSize    = 18;
    H = ylabel('RN [e]');
    H.Interpreter = 'latex';
    H.FontSize    = 18;

    %print LimMagSN10_R3000_sampling0.8.jpg -djpeg100
    
end




