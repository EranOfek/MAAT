function [D,S,Sc,PDn,PDr]=imsub_core(SimN,SimR,varargin)
%--------------------------------------------------------------------------
% imsub_core function                                           class/@SIM
% Description: The core functionality of proper image subtraction.
%              Given a SIM image, a reference image and their PSFs,
%              calculate the proper image subtraction.
%              This function assumes that:
%              1. The images gain are 1.
%              2. The image include the background.
%              3. The background and std are available in the 'BackIm' and
%                 'ErrIm' fields of the SIM, respectively.
%              4. The PSF is available in the 'PSF' field of the SIM.

%              If needed, flux matching maybe done by this function.
% Input  : - SIM images of the new.
%            The image should be background subtracted and registered to
%            the reference image.
%            The original StD and Background images may be provided in
%            the 'ErrIm' and 'BackIm' fields.
%          - A SIM image containig the reference image.
%            The image should be background subtracted.
%            The original StD and Background images may be provided in
%            the 'ErrIm' and 'BackIm' fields.

%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:

% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

ImageField       = 'Im';
BackField        = 'BackIm';
ErrField         = 'ErrIm';
PSFField         = 'PSF';



DefV.ImageField         = ImageField;
DefV.BackField          = BackField;
DefV.ErrField           = ErrField;
DefV.F_r                = 1;
DefV.Epsilon            = 1e-8;
DefV.CalcDback          = true;
DefV.DbackPar           = {'Block','full'};
DefV.CalcSback          = true;
DefV.SbackPar           = {'Block','full'};
DefV.MeanStdFun         = @mean;
DefV.FftShift           = true;
DefV.ReadParPSF         = {};
DefV.SigmaPosX          = [];
DefV.SigmaPosY          = [];
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


ImageField = InPar.ImageField;
BackField  = InPar.BackField;
ErrField  = InPar.ErrField;


NsimN = numel(SimN);
NsimR = numel(SimR);
Nsim  = max(NsimN,NsimR);

InPar.F_r = InPar.F_r(:).*ones(NsimR,1);

if ~(NsimR==1 || NsimR==NsimN),
    error('Reference SIM must have number of elements equal 1 or that of New SIM');
end


% SigmaPosX/Y can be:
% A scalar
% A vector (element per image)
% A SIM
if (SIM.issim(InPar.SigmaPosX))
    SigmaPosX = sim2cell(InPar.SigmaPosX);
else
    SigmaPosX = num2cell(InPar.SigmaPosX);
end
if (SIM.issim(InPar.SigmaPosY))
    SigmaPosY = sim2cell(InPar.SigmaPosY);
else
    SigmaPosY = num2cell(InPar.SigmaPoY);
end
Nsigposx = numel(SigmaPosX);
Nsigposy = numel(SigmaPosY);



CalcD = true;
if (nargin>1),
    CalcS = true;
    if (nargin>2),
        CalcSc = true;
    end
end

% get the PSFs
NewPSF = getmpsf(SimN);
RefPSF = getmpsf(SimR);


% if (CalcD),
D   = SIM(Nsim,1);
PDr = ClassPSF(Nsim,1);
PDn = ClassPSF(Nsim,1);

% end
if (CalcS),
    S = SIM(Nsim,1);
end
if (CalcSc),
    Sc = SIM(Nsim,1);
end
for Isim=1:1:Nsim,
    %----------------------------
    %--- For each SIM element ---
    %----------------------------
    
    % Indices for SimN and SimR:
    IsimN = min(Isim,NsimN);
    IsimR = min(Isim,NsimR);
    
    %------------------------------------
    %--- Read the images from the SIM ---
    %------------------------------------
    % Get images
    NewWB    = SimN(IsimN).(ImageField);    % New Image (with background)
    RefWB    = SimR(IsimR).(ImageField);    % Ref Image (with background)
    % Get background
    NewBack  = SimN(IsimN).(BackField);     % New background
    RefBack  = SimR(IsimR).(BackField);     % Ref background
    % Get std 
    NewStd   = SimN(IsimN).(ErrField);      % New back std
    RefStd   = SimR(IsimR).(ErrField);      % Ref back std
    % Background subtracted images
    NewBS    = NewWB - NewBack;             % Background subtracted New image
    RefBS    = RefWB - RefBack;             % Background subtracted Ref image
  
    % Image size
    NewSize = size(NewWB);
    RefSize = size(RefWB);
    if (~all(NewSize==NewRef)),
        error('New image number %d and reference image number %d do not have the same size',IsimN,IsimR);
    end

    
    % Read the PSF and convert it to have the same size as the image
    % including shift to the origin
    PSF_N = psf_pad_fftshift(NewPSF{IsimN},NewSize,InPar.FftShift,InPar.ReadParPSF{:});
    PSN_R = psf_pad_fftshift(RefPSF{IsimR},RefSize,InPar.FftShift,InPar.ReadParPSF{:});
    
    % renormalize PSF (just in case)
    PSF_N    = PSF_N./sum(PSF_N(:));
    PSF_R    = PSF_R./sum(PSF_R(:));

    % fft2 of images
    fftPn   = fft2(PSF_N);    % FFT of P_n
    fftPr   = fft2(PSF_R);    % FFT of P_r
    fftN    = fft2(NewBS);    % FFT of N (back. subtracted)
    fftR    = fft2(RefBS);    % FFT of R (back. subtracted)
    fftNwb  = fft2(NewWB);    % FFT of N (not back. sub.)
    fftRwb  = fft2(RefWB);    % FFT of R (not back. sub.)
    
    FPrFN = fftPr.*fftN;   % fft(P_r).*fft(N)
    FPnFR = fftPn.*fftR;   % fft(P_n).*fft(R)

    if (isempty(InPar.Beta)),
        %---------------------
        %--- Estimate Beta ---
        %---------------------
       
    else
        % Beta provided by user
        Beta = InPar.Beta;
    end
    
    %------------------------
    %--- Flux definitions ---
    %------------------------
    % F_r is the flux-based zero point of the Ref image.
    % default is 1.
    F_r = InPar.F_r(IsimR);  
    % F_n is the flux-based zero point of the New image.
    F_n = F_r.*Beta;
    
    %-------------------------
    %--- Noise definitions ---
    %-------------------------
    sigmaN = InPar.MeanStdFun(NewStd);
    sigmaR = InPar.MeamStdFun(RefStd);
    
    %--------------------------------------------------
    %--- Calculate the proper subtraction image (D) ---
    %--------------------------------------------------
    FPn2 = FPn.*conj(FPn);       % |fft(P_n)|^2
    FPr2 = FPr.*conj(FPr);       % |fft(P_r)|^2
    % The denominator of S (Equation 12 in ZOGY 2016)   
    % Epsilon is added in order to avoid machin precission problems
    % Doesn't influence the results
    % Don't use 1/sqrt() as this will introduce division by zero...
    S_denom = sigmaR.^2.*F_n.^2.*FPn2 + sigmaN.^2.*F_r.^2.*FPr2 + InPar.Epsilon;
    D_denom = sqrt(S_denom);
    % Equation 13 in ZOGY 2016:
    fftD    = (F_r.*fftPr.*fftN - F_n.*fftPn.*fftR)./D_denom;
    F_denom = sqrt(sigmaN.^2.*F_r.^2 + sigmaR.^2.*F_n.^2); 
    F_D     = F_r.*F_n./F_denom;
    fftP_D  = F_r.*F_n.*fftPr.*fftPn./(F_D.*D_denom);
    
    % Claculate P_Dn and P_Dr
    % Equations 18, 19, E13, E16, E17 in ZOGY 2016
    F_Dn    = F_r./F_denom;
    F_Dr    = F_n./F_denom;
    fftP_Dn = F_r.*fftPr./(F_Dn.*D_denom);
    fftP_Dr = F_n.*fftPn./(F_Dr.*D_denom);

    
    %--------------------------------------------
    %--- Populate the proper image (D) output ---
    %--------------------------------------------
    D(Isim).(ImageField)  = fftshift(ifft2(fftD));
    D(Isim).(PSFField)    = fftshift(ifft2(fftP_D));
    PDn(Isim).(PSFField)  = fftshift(ifft2(fftP_Dn));
    PDr(Isim).(PSFField)  = fftshift(ifft2(fftP_Dr));

    if (InPar.CalcDback),
        % calculate the back and std of D
        D(Isim) = background(D(Isim),InPar.DbackPar{:});
    end
    
    %---------------------------------------------------------------------
    %--- Calculate the Score image statistics (S) for source detection ---
    %---------------------------------------------------------------------
    % S is given by filtering D with P_D
    if (CalcS),
        fftS    = F_D.*fftD.*conj(fftP_D);

        %-------------------------------------------
        %--- Populate the score image (S) output ---
        %-------------------------------------------
        S(Isim).(ImageField)  = fftshift(ifft2(fftS));
        S(Isim).(PSFField)    = D(Isim).(PSFField);
        if (InPar.CalcSback),
            % calculate the back and std of S
            S(Isim) = background(S(Isim),InPar.SbackPar{:});
        end
        
        %---------------------------
        % Transient flux estimation
        %---------------------------
        % Flux normalization (F_S):
        % the abs value is needed in order to get rid of small imaginary part
        % Equation 42 in ZOGY 2016
        F_S = abs(sumnd((F_n.^2.*FPn2 .* F_r.^2.*FPr2)./S_denom));
        
       
    end
    
    %-------------------------------------------------------------------------
    %--- Calculate the corrected score image (S_corr) for source detection ---
    %-------------------------------------------------------------------------
    if (CalcSc),
        
            
        if (isempty(InPar.SigmaPosX) || isempty(InPar.SigmaPosY)),
            %Vadd = 0;
            FgradN = 0;
            FgradR = 0;
        else
            IsigposX = min(IsimN, Nsigposx);
            IsigposY = min(IsimN, Nsigposy);

            Kn  = conj(fftPn).*FPr2./S_denom;  % This is not the Kn definition in the paper (missing Fn,Fr)
            Kr  = conj(fftPr).*FPn2./S_denom;  % This is not the Kr definition in the paper (missing Fn,Fr)

            % Kn is in Fourier space - the gradient is calculated in real space
            [GradNx,GradNy] = gradient(ifft2(F_n         .*Kn.*fftN));  % gradient of S_n
            [GradRx,GradRy] = gradient(ifft2(F_r.*Beta.^2.*Kr.*fftR));  % gradient of S_r
            %[GradNx,GradNy] = gradient(ifft2(Kn.*FN));  % gradient of S_n      <-----------------
            %[GradRx,GradRy] = gradient(ifft2(Kr.*FR));  % gradient of S_r      <-----------------
            FgradN = (SigmaPosX{IsigposX}.*GradNx).^2 + (SigmaPosY{IsigposY}.*GradNy).^2;
            FgradR = (SigmaPosX{IsigposX}.*GradRx).^2 + (SigmaPosY{IsigposY}.*GradRy).^2;

            %VaddFT = Kn2.*FgradN + Kr2.*FgradR;
            % Any additional noise (e.g., color-refraction noise) should be
            % added to Vadd.
            %Vadd = FgradN + FgradR;

        end

        % Variance correction
        % Incluing: Source noise and astrometric noise
        %Sub(Iim).Snoise = sqrt(ifft2(Kn2.*FNnbs + Kr2.*FRnbs) + Vadd);
        Snoise = sqrt(ifft2(Kn2.*FNnbs + Kr2.*FRnbs) + FgradN + FgradR);

        
        % estimate the ratio between N and R noise:
        %sqrt(mediannd((ifft2(Kn2.*FNnbs) + FgradN)./(ifft2(Kr2.*FRnbs) + FgradR)))

        % Error in F_S flux measurment parameter
        % Equation 43 in ZOGY 2016
        Std_S = Snoise./F_S;

    end
end

    



