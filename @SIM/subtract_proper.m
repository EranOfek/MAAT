function [Summary,D,S,Scorr,SigmaF,GradSN]=subtract_proper(Sim,SimRef,varargin)
% Image subtraction using the ZOGY algorithm.
% Package: class/@SIM
% Description: Image subtraction using the Zackay, Ofek, Gal-Yam (2016)
%              method. 
% Input  : - A SIM object. If the Sim object does not contain a background,
%            error, and PSF then they will be generated.
%          - A single element SIM object containing a reference image.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'GainCorrect' - A flag indicating if to multiply the input
%                       images by the gain factor. Default is false.
%                       Note that in order for ZOGY to work properly the
%                       images needed to be multiplied by
%                       their gain (i.e., having electron units)
%            'GainCorrectRef' - A flag indicating if to multiply the
%                       reference image by the gain factor.
%                       Default is false.
%            'BackReCalc' - A flag indicating if to recalculate the
%                       background and error.
%                       Default is true.
%            'BackPar' - A cell array of additional parameters to pass to
%                       the background function.
%                       Default is {'Block','full'}.
%                       Note that this version of ZOGY requires a global
%                       background and error estimate per image.
%            'FluxMatchMethod' - Flux match method:
%                       'zp' : calculate the flux matching using PSF
%                              photometry of the new and reference images.
%                       'robustfit' - The ZOGY original method using
%                              robust fitting.
%                       'ransac' - The ZOGY original method using
%                              ransac fitting.
%                       'known' - value is provided by user.
%                       Default is 'robustfit'.           
%            'FluxMatch' - A vector of user provided flux matching value
%                       per image. If given then 'FluxMatchMethod' will
%                       be set to 'known'.
%                       Default is [].
%            'BetaTol' - Tolerance for flux matching convergance using
%                       the robustfit and ransac methods.
%                       Default is 1e-4.
%            'RelZPPar'- Cell array of additional paramaeters to pass to
%                       rel_zp.m for calculating the relative zero point
%                       (flux mataching). Default is {}.
%            'PN'     - A matrix of ClassPSF object containing the PSF for
%                       the new image. If given this will override the PSF
%                       provided in the SIM object.
%                       Default is [].
%            'PR'     - The same as 'PN', but for the reference image.
%                       Default is [].
% 
% DefV.PopPSFPar          = {};
% DefV.Align              = false;
% DefV.AlignPar           = {};
% DefV.AlignFields        = {SIM.ImageField, SIM.BackField, SIM.ErrField, MASK.MaskField};
% DefV.FunGlobalErr       = @median;   % Fun defined on vector
% DefV.Eps                = 1e-8;      % constant for avoiding division by zero - not really necessery
% 
% DefV.RF_Wfun            = 'welsch';  %{}; %'fair';     % robust fit Wfun parameter - set to {} in order to use default.
% DefV.MaxIterBeta        = 10;        % maximum number of iterations for beta solver
% DefV.UpdateGamma        = false;      % subtrcat updated background from new
% DefV.SigmaX             = [];        % X astrometric error [pix] % scalar, vector, SIM
% DefV.SigmaY             = [];        % Y astrometric error [pix] % scalar, vector, SIM
% DefV.SigmaField         = SIM.ImageField;
% 
% DefV.ReNormScorr        = true;
% DefV.ReNormBlock        = [256 256];

% DefV.Verbose            = true;
% DefV.VerbosePlot        = true;

%
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Su,D,S,Scorr]=subtract_proper(Sim,SimRef);
%          [Su,D,S,Scorr]=subtract_proper(N,R,'GainCorrect',false,'GainCorrectRef',false,'FluxMatchMethod','robustfit','SigmaX',0.02,'SigmaY',0.02,'Align',false)
% Reliable: 2
%--------------------------------------------------------------------------

ImageField   = SIM.ImageField;   % 'Im';
BackField    = SIM.BackField;
ErrField     = SIM.ErrField;     % 'ErrIm';
WeightField  = SIM.WeightField;  % 'WeightIm';
MaskField    = SIM.MaskField;    % 'Mask';
PSFField     = 'PSF';
%MaskDicField = 'MaskDic';

DefV.GainCorrect        = false;
DefV.GainCorrectRef     = false;
DefV.BackReCalc         = true;
DefV.BackPar            = {'Block','full'};
DefV.FluxMatchMethod    = 'robustfit';    % 'zp'|'robustfit'|'ransac'|'known'
DefV.FluxMatch          = [];      % If provided than override over FluxMatchMethod
DefV.BetaTol            = 1e-4;
DefV.RelZPPar           = {};

DefV.PN                 = [];      % matrix or ClassPSF with PSF to override
DefV.PR                 = [];

DefV.PopPSFPar          = {};
DefV.Align              = false;
DefV.AlignPar           = {};
DefV.AlignFields        = {SIM.ImageField, SIM.BackField, SIM.ErrField, MASK.MaskField};
DefV.FunGlobalErr       = @median;   % Fun defined on vector
DefV.Eps                = 1e-8;      % constant for avoiding division by zero - not really necessery

DefV.RF_Wfun            = 'welsch';  %{}; %'fair';     % robust fit Wfun parameter - set to {} in order to use default.
DefV.MaxIterBeta        = 10;        % maximum number of iterations for beta solver
DefV.UpdateGamma        = false;      % subtrcat updated background from new
DefV.SigmaX             = [];        % X astrometric error [pix] % scalar, vector, SIM
DefV.SigmaY             = [];        % Y astrometric error [pix] % scalar, vector, SIM
DefV.SigmaField         = SIM.ImageField;

DefV.ReNormScorr        = true;
DefV.ReNormBlock        = [256 256];



DefV.Verbose            = true;
DefV.VerbosePlot        = true;

% 
% 
% DefV.ExecField          = SIM.ImageField;
% DefV.GainCorrect        = true;
% DefV.BackReCalc         = false;
% DefV.BackPar            = {'Block',[256 256]};
% DefV.FunGlobalBack      = @median;   % Fun defined on SIM and of the form: median(Sim,'ExecField',ErrField)
% DefV.WeightVar          = [];
% DefV.WeightTran         = [];
% DefV.PhotZPPar          = {};
% DefV.FilterPar          = {};
% DefV.CoaddPar           = {};

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (nargout>3)
    CalcScorr = true;
else
    CalcScorr = false;
end
if (nargout>2 || CalcScorr)
    CalcS = true;
else
    CalcS = false;
end



Nsim  = numel(Sim);
NsigmaX = numel(InPar.SigmaX);
NsigmaY = numel(InPar.SigmaY);

%------------------------------
%--- Populate the meta data ---
%------------------------------
% populate the images with all the necessery information:
% Gain correct
if (InPar.GainCorrect)
    Sim    = gain_correct(Sim);
end
if (InPar.GainCorrectRef)
    SimRef = gain_correct(SimRef);
end

% Background and error image
if (InPar.BackReCalc)
    BackExist    = false(Nsim,1);
    BackExistRef = false(1,1);
else
    BackExist    = isfield_populated(Sim,BackField);
    BackExistRef = isfield_populated(SimRef,BackField);
end
Sim(~BackExist)       = background(Sim(~BackExist),      InPar.BackPar{:});
SimRef(~BackExistRef) = background(SimRef(~BackExistRef),InPar.BackPar{:});

% Catalog [if transperancy using catalog is needed]
if (~isempty(InPar.FluxMatch))
    % FluxMatch is given
    InPar.FluxMatchMethod = 'known';
else
    % FluxMatch is unknown
    switch lower(InPar.FluxMatchMethod)
        case 'zp'
            % estimate flux-match (transperancy) via relative photometry
            SimAll = [SimRef; Sim(:)];  % All SIM including ref image
            [ResZP,TranZP,SimAll] = phot_zp(SimAll,InPar.RelZPPar{:},'MatchPar',{'SkipWCS',true});
            % update the flux based transperancies
            % FluxMatch is calculated relative to reference image
            InPar.FluxMatch = TranZP(2:end)./TranZP(1);
            Sim    = SimAll(2:end);
            SimRef = SimAll(1);
            if (InPar.Verbose)
                fprintf('Flux match ratio: %f\n',InPar.FluxMatch);
            end
            % Set FluxMatchMethod to known
            InPar.FluxMatchMethod = 'known';
        otherwise
            % do nothing
    end
end

% Populate the PSF
Sim    = populate_psf(Sim,InPar.PopPSFPar{:});
SimRef = populate_psf(SimRef,InPar.PopPSFPar{:});

%------------------------------------------------
%--- Align images relative to reference image ---
%------------------------------------------------
if (InPar.Align)
    [AlSim,ResFit,IndRef] = align(Sim,SimRef,InPar.AlignPar{:},'ExecField',InPar.AlignFields);
end


%-------------------------------------
%--- Perform the image subtraction ---
%-------------------------------------
% Set the reference image and its data
Rwb    = SimRef.(ImageField);   % R with background
Rback  = SimRef.(BackField);
Rsigma = InPar.FunGlobalErr(SimRef.(ErrField)(:));
R      = Rwb - Rback;

% check that R doesn't contain NaNs
if (any(isnan(R)))
    error('Reference image contains NaNs');
end

Rsize  = fliplr(imagesize(SimRef));  % [Y X] imagesize
NpixTotal = prod(Rsize);
% PSF, padding and shift to origin
PR     = psf_pad2mat(SimRef,Rsize,true);
% fft of R and PR
fftR   = fft2(R);
fftPR  = fft2(PR);

D      = SIM(size(Sim));
S      = SIM(size(Sim));
Scorr  = SIM(size(Sim));
SigmaF = SIM(size(Sim));
%VSc    = SIM(size(Sim));

Nsim = numel(Sim);
for Isim=1:1:Nsim
    % For each image in SIM
    
    % set the new image and its data
    Nwb    = Sim(Isim).(ImageField);   % N with background
    % make sure Nwb do not contain negative values
    % NO IDEA WHY I ADDED THIS LINE!
    %Nwb(Nwb<0) = 0;
    
    Nback  = Sim(Isim).(BackField);
    Nsigma = InPar.FunGlobalErr(Sim(Isim).(ErrField)(:));
    N      = Nwb - Nback;
    if (any(isnan(N)))
        error('New image %d contains NaNs',Isim);
    end

    Nsize  = fliplr(imagesize(Sim(Isim)));  % [Y X] imagesize
    % PSF, padding and shift to origin
    PN     = psf_pad2mat(Sim(Isim),Nsize,true);
    % fft of N and PN
    fftN   = fft2(N);
    fftPN  = fft2(PN);

    
    % Image subtraction factors
    % This factors are indepndent of the flux matching (with the
    % exception of N that may change if gamma may change):
    fftPRfftN = fftPR.*fftN;   % fft(P_r).*fft(N)
    fftPNfftR = fftPN.*fftR;   % fft(P_n).*fft(R)

    %------------------
    %--- Flux Match ---
    %------------------
    % Check if FluxMatch is still needed
    switch lower(InPar.FluxMatchMethod)
        case 'known'
            % Beta is the FluxMatch factor for the current image
            Beta     = InPar.FluxMatch(Isim);
            Gamma    = 0;
            IterBeta = NaN;
        case 'min_scorr'
            % Flux matching by minimizing the residuals in Scorr
            Nf=filter(sub_background(Sim),[],'NormErr',true);
            FlagSources = abs(Nf.(ImageField))>10 & Sim(Isim).(ImageField)<50000 & SimRef.(ImageField)<50000;
            
            Beta0   = 1;
            BetaVec = [0.98, 1.00, 1.05]; 
            %BetaVec = [0.9:0.01:1.1];
            Nbv     = numel(BetaVec);
            
            MaxLoop = 5;
            Loop = 0;
            Converge = false;
            MaxLoop = 3;
            while (~Converge && Loop<MaxLoop)
                Loop = Loop + 1;
                
                for Ibv=1:1:Nbv
                    [~,~,~,Scorr] = subtract_proper(Sim(Isim),SimRef,...
                                                    'GainCorrect',false,'GainCorrectRef',true,...
                                                    'SigmaX',InPar.SigmaX(Isim),'SigmaY',InPar.SigmaY(Isim),...
                                                    'Align',false,...
                                                    'FluxMatch',BetaVec(Ibv).*Beta0);
                    BetaStd(Ibv) = Util.stat.rstd(Scorr(FlagSources));
                    
                    
                    %BetaStd(Ibv) = mean(Scorr(FlagSources));
                end

                plot(BetaVec,BetaStd,'o')
                
                P = polyfit(BetaVec,BetaStd,2);
                BetaPrev = Beta0;
                if (P(1)>0)
                    % minimum found
                    Beta0 = -P(2)./(2.*P(1));
                else
                    Beta0 = 2.*Beta0 + P(2)./(2.*P(1));
                end
                if (abs(Beta0-BetaPrev)<0.03)
                    Converge = true;
                end
                
            end
            Beta = Beta0
            

            
            
        case {'robustfit','ransac'}
            % Use the FluxMatch method suggested in the original
            % ZOGY paper -- find beta and gamma
                       
            % solve for beta and gamma (no translations)
            % using iterative robust fitting

            % set initial guess for the value of beta and gamma
            Beta     = 1;
            GammaTag = 0;

            BetaGammaConverg = false;
            IterBeta = 0;
            while (~BetaGammaConverg && IterBeta<=InPar.MaxIterBeta)
                % The denominator of D with the beta term:
                IterBeta = IterBeta + 1;
                Denominator     = Beta.^2.*Rsigma.^2 .* fftPN.*conj(fftPN) + Nsigma.^2 .* fftPR.*conj(fftPR) + InPar.Eps;
                SqrtDen = sqrt(Denominator);
                % ZOGY: Equation 37
                Dn = ifft2(fftPRfftN./SqrtDen);
                % ZOGY: Equation 38
                Dr = ifft2(fftPNfftR./SqrtDen);

                % select only pixels which are >1 sigma above the image
                % noise.
                % Is this necesssery?
                [Back_Dn,Std_Dn] = Util.stat.mode_fit(Dn);
                [Back_Dr,Std_Dr] = Util.stat.mode_fit(Dr);
                %Il0 = Dn(:)>1 & Dr(:)>1 & InPar.SatPixR(:)==0 & InPar.SatPixN(:)==0;
                %Il0 = Dn(:)>Std_Dn.*7 & Dr(:)>Std_Dr.*7 & InPar.SatPixR(:)==0 & InPar.SatPixN(:)==0;

                % FFU - remove also bad/saturated pixels according
                % to mask
                Nsig = 10;
                Il0 = Dn(:)>Std_Dn.*Nsig & Dr(:)>Std_Dr.*Nsig & Nwb(:)<40000; %& Rwb(:)<50000;
                %Il0 = Dn(:)>InPar.SigmaN(IsigN) & Dr(:)>InPar.SigmaR & InPar.SatPixR(:)==0 & InPar.SatPixN(:)==0;
                %Il0 = Dn(:)>-Inf & Dr(:)>-Inf;

                % robust fitting of the Dn and Dr

                switch lower(InPar.FluxMatchMethod)
                    case 'robustfit'
                        Par = robustfit(Dr(Il0),Dn(Il0),InPar.RF_Wfun);
                    case 'ransac'
                        Par = Util.fit.fitlin_ransac(Dr(Il0),Dn(Il0),Dn(Il0),1,100);
                    otherwise
                        error('Unknown BetaMethod option');
                end
                GammaTag = Par(1);
                OldBeta  = Beta;
                Beta     = Par(2);

                BetaGammaConverg = abs(OldBeta - Beta)<InPar.BetaTol;

            end
            Gamma = GammaTag.*sqrt(Nsigma.^2 + Beta.^2.*Rsigma.^2);
            % update the background level of the New
            if (InPar.UpdateGamma)
                % update N according to best fit gamma shift
                N         = N - Gamma;
                fftN      = fft2(N);
                fftPRfftN = fftPR.*fftN;
            end

            if (InPar.VerbosePlot)
                % plot Dr vs. Dn and best fit line
                plot(Dr(Il0),Dn(Il0),'.');
                hold on
                plot(Dr(Il0),Par(1) + Par(2).*Dr(Il0),'r.')
            end

            if (InPar.Verbose)
                fprintf('Flux matching using %s method\n',InPar.FluxMatchMethod);
                fprintf('   Solution found after %d iterations\n',IterBeta);
                fprintf('   gamma_tag (offset) = %f\n',Par(1));
                fprintf('   beta (N/R gain)    = %f\n',Par(2));
            end

              
        otherwise
            error('Unknown FluxMatchMethod option');
    end  % end of FluxMatch switch
    
    % continue image subtraction
    % Define the flux ratio for image R and N
    FR = 1;      % set the flux of the reference to 1
    FN = Beta;   % Beta=F_N/F_R - ZOGY Equation 35
    
    Summary(Isim).Beta = Beta;
    
    % |fft(P_n)|^2
    fftPN2 = fftPN.*conj(fftPN);       
    % |fft(P_r)|^2
    fftPR2 = fftPR.*conj(fftPR);       
    % Dnominator of ZOGY Equation 12
    Denominator     = Nsigma.^2 .* FR.^2 .* fftPR2 + Rsigma.^2 .* FN.^2 .* fftPN2 + InPar.Eps;
    SqrtDen         = sqrt(Denominator);             % don't use 1/sqrt() as this will introduce division by zero...
    
    %-------------------------------------------------------------------
    %--- Calculate the de-corelated image subtraction statistics (D) ---
    %-------------------------------------------------------------------
    % ZOGY Equation 13
    %------------------
    % The Image subtraction statistics D
    % This is the noise decorlated image subtraction
    fftD = (FR.*fftPRfftN - FN.*fftPNfftR)./SqrtDen;
    D(Isim).(ImageField) = ifft2(fftD);
    
    % ZOGY Equation 15
    FD = FR.*FN./sqrt(Nsigma.^2 .* FR.^2 + Rsigma.^2 .* FN.^2);
    Summary(Isim).FD = FD;
    
    % ZOGY Equation 14
    %------------------
    % The PSF P_D of the difference statistics D:
    fftPD = FR.*FN .* fftPR.*fftPN./(FD .* SqrtDen);
    
    D(Isim).(PSFField) = fftshift(ifft2(fftPD));
    
    % ZOGY Equation 18
    %  FFU
    % ZOGY Equation 19
    %  FFU
    
    
    
    %---------------------------------------------------------------------
    %--- Calculate the Score image statistics (S) for source detection ---
    %---------------------------------------------------------------------
    if (CalcS)
        % ZOGY Equation 28
        fftKR  = FR .* FN.^2 .* conj(fftPR).*fftPN2./Denominator;
        % ZOGY Equation 29
        fftKN  = FN .* FR.^2 .* conj(fftPN).*fftPR2./Denominator;          % This is not the Kr definition in the paper (missing Fn,Fr)

        % ZOGY Equation 12
        %------------------
        % The filtered image subtraction statistics S:
        % This is the log-likelihood to find a transient (with the image PSF)
        % at each pixel.
        % alternatively can calculated by filtering D with P_d:    
        S(Isim).(ImageField) = ifft2(fftKN.*fftN - fftKR.*fftR);
        
        
        % Alternatively (the same answer)
        % ZOGY Equation 17
        %S(Isim).(ImageField) = FD.*ifft2(fftD.*conj(fftPD));
    end
    
    %-----------------------
    %--- Calculate Scorr ---
    %-----------------------
    % Calculate the source noise and astrometric noise corrected
    % transient detection statistics (Scorr): ZOGY Equation 25
    if (CalcScorr)
        % Note that Kr and Kn are calculated in the "CalcS" block:
        % fftKR2 is fft(k_r^2) needed in ZOGY Equation 27:
        fftKR2 = fft2(ifft2(fftKR).^2);   
        % fftKN2 is fft(k_n^2) needed in ZOGY Equation 26:
        fftKN2 = fft2(ifft2(fftKN).^2);   
        
                
        
        if (isempty(InPar.SigmaX) || isempty(InPar.SigmaY))
            VastSN = 0;
            VastSR = 0;
            GradSNx = 0;
            GradSNy = 0;
            GradSRx = 0;
            GradSRy = 0;
            VastSN  = 0;
            VastSR  = 0;
        else
            IndX     = min(Isim,NsigmaX);
            IndY     = min(Isim,NsigmaY);
            if (SIM.issim(InPar.SigmaX(IndX)))
                % InPar.SigmaX is a SIM image
                % SigmaX is an image:
                SigmaX = InPar.SigmaX(IndX).(InPar.SigmaField);
            else
                % SigmaX is a scalar:
                SigmaX = InPar.SigmaX(IndX);
            end
            if (SIM.issim(InPar.SigmaY(IndY)))
                % InPar.SigmaY is a SIM image
                % SigmaY is an image:
                SigmaY = InPar.SigmaY(IndY).(InPar.SigmaField);
            else
                % SigmaY is a scalar:
                SigmaY = InPar.SigmaY(IndY);
            end
            
            % calculate the gradient of S_N
            % needed for ZOGY Equation 30
            % S_N is defined in ZOGY Equation 31
            [GradSNx,GradSNy] = gradient(ifft2(fftKN.*fftN));  % gradient of S_n
            % S_R is defined in ZOGY Equation 33
            [GradSRx,GradSRy] = gradient(ifft2(fftKR.*fftR));  % gradient of S_n
            
            % ZOGY Equation 30
            VastSN = (SigmaX.*GradSNx).^2 + (SigmaY.*GradSNy).^2;

            % ZOGY Equation 32
            VastSR = (SigmaX.*GradSRx).^2 + (SigmaY.*GradSRy).^2;
            
        end
        % Any additional noise (e.g., color-refraction noise) should be
        % added to Vadd here after the VastSN, VastSR terms
        % Denominator of ZOGY Equation 25
        VScorr = sqrt(ifft2(fftKN2.*fft2(Nwb) + fftKR2.*fft2(Rwb)) + VastSN + VastSR);
        
        %[MeanS,StdS] = Util.stat.mode_fit(S);
        %S = S - MeanS;
        Scorr(Isim).(ImageField) = S(Isim).(ImageField)./VScorr;
        
        %VSc(Isim).(ImageField) = VScorr;
        
        %-------------------------
        %--- renormalize Scorr ---
        %-------------------------
        % If the images have gradients in their background and or noise,
        % Scorr will have gradients and its local mean will deviate from zero,
        % and local std will deviate from one.
        % This section calculate the local background and noise of Scorr
        % and renormalize it
        % Note that this is an approximation to solving D locally (i.e.,
        % pixel by pixel).
        if (InPar.ReNormScorr)
            % Renormalize Scorr
            if any(~isreal(Scorr(Isim).(ImageField)))
                warning(sprintf('Scorr of image %d contains imaginary numbers',Isim));
                Scorr(Isim) = abs(Scorr(Isim));
            end
            Scorr(Isim) = background(Scorr(Isim),'Block',InPar.ReNormBlock);
            Scorr(Isim).(ImageField) = (Scorr(Isim).(ImageField) - Scorr(Isim).(BackField)); %./Scorr(Isim).(ErrField);
            RStd = Util.stat.rstd(Scorr(Isim).(ImageField)(:));
            Scorr(Isim).(ImageField) = Scorr(Isim).(ImageField)./RStd;
            Scorr(Isim).(ErrField)   = RStd;
        end
        
        
        %---------------------------
        % Transient flux estimation
        %---------------------------
        % ZOGY Equation 42
        % Flux normalization (F_S):
        % the abs value is needed in order to get rid of small imaginary part
        % Note that the division by NpixTotal is needed because
        % fft(ones)=sum(ones)
        F_S = (Util.stat.sumnd((FN.^2 .* fftPN2 .* FR.^2 .* fftPR2)./Denominator))./NpixTotal;
        if (imag(F_S)>0)
            warning('Negative F_S for image %d - use abs value',Isim);
            F_S = abs(F_S);
            
        end
        
        %F_S = abs(Util.stat.sumnd(FN.^2.*FR.*fftPN2.*conj(fftPR)./Denominator))
%         AAA = FN.^2.*FR.*fftPN2.*conj(fftPR)./Denominator;
%         F_S = AAA(1,1);
%         
        % Flux error : ZOGY Equation 43
        % assumes images are in units of electrons or gain is known
        SigmaFlux = VScorr./F_S;

        SigmaF(Isim).(ImageField) = SigmaFlux;
        
        Summary(Isim).F_S = F_S;
        Summary(Isim).SigmaFlux = SigmaFlux;
        
    end
    GradSN(Isim).X = GradSNx;
    GradSN(Isim).Y = GradSNy;
    
    
    
    
end
