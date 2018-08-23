function [Sub]=imsub_fft(Ref,New,varargin)
% Proper image subtraction - OBSOLETE: user SIM/subtract_proper instead.
% Package: ImUtil.Im
% Description: Perform optimal proper image subtraction using the
%              method of Zackay, Ofek & Gal-Yam (2016).
%              Given a reference image and a new image, return the
%              proper image subtraction D, the uncorrected score S,
%              the score correction, the corrected score,
%              the PSF of D, the PSF of a delta
%              function in N as appear in D, and the PSF of a delta
%              function in R as appear in D.
%              The program can take care of the image registration,
%              background subtraction, variance estimation,
%              and PSF measurment.
%              OBSOLETE: user SIM/subtract_proper instead.
% Input  : - A single reference image. This can be a FITS file, a SIM
%            image or any format supported by image2sim.m.
%          - A list of new images. This can be FITS files, a SIM
%            images or any format supported by images2sim.m.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'GainN'  - Gain [e-/ADU] of new image. This is either the
%                       gain value or a string containing the gain
%                       FITS header keyword. Default is 'GAIN'.
%            'GainR'  - Gain [e-/ADU] of ref image. This is either the
%                       gain value or a string containing the gain
%                       FITS header keyword. Default is 'GAIN'.
%            'Align'  - Align images prior to subtraction {true|false}.
%                       Default is false.
%            'SubBackN'-Subtract background from new image {true|false}.
%                       Default is true.
%            'SubBackR'-Subtract background from ref image {true|false}.
%                       Default is true.
%            'BackMethod'-Background subtraction method {'mode_fit'}.
%                       Default is 'mode_fit' (see mode_fit.m for details).
%            'SigmaN' - Sigma noise of new image. If empty, then will
%                       estimate sigma using mode_fit.m.
%                       Default is [].
%            'SigmaR' - Sigma noise of ref image. If empty, then will
%                       estimate sigma using mode_fit.m.
%                       Default is [].
%            'MethodStD'- Sigma estimation method.
%                      If 'fit' then use mode_fit.m.
%                      Default is 'fit'.
%            'RN_N'   - Readnoise [e-] of new image. Default is 5.
%            'RN_R'   - Effective readnoise [e-] of ref image. Default is 5.
%            'PSF_N'  - A matrix, or a cell arrya of matrices (per each new
%                       image) of PSFs in the new image. If empty then
%                       will attempt to estimate the PSF.
%                       The PSF should be centered, otherwise, this will
%                       introduce a shift to the subtraction image.
%                       Default is empty.
%            'PSF_R'  - A matrix of the PSF in the ref image. If empty then
%                       will attempt to estimate the PSF.
%                       The PSF should be centered, otherwise, this will
%                       introduce a shift to the subtraction image.
%                       Default is empty.
%            'PSFbuildPar' - Cell array of key,val parameters to pass to
%                       psf_builder.m. Default is
%                       {'BoxHalfSize',[15
%                       15],'Rad0',15,'ZeroEdgesThresh',2,'ZeroEdgesMethod','rad'}.
%            'Beta'     - Value of  beta (zero point matching between the
%                       images). If empty then will attempt to estimate beta.
%                       Default is empty.
%            'BetaMethod' - Method to use for estimating beta (zero point
%                       matching between the images). Options are:
%                       'robustfit' - use robust fitting of Dn and Dr.
%                                     Default.
%            'BetaTol' - Tolerance for estimating beta.
%                       Default is 1e-4.
%            'RF_Wfun'- robustfit Wfun parameter. Set to {} in order to 
%                       use default. Default is {}.
%            'MaxIterBeta'- Maximum number of iterations for beta solver.
%                       Default is 10.
%            'UpdateGamma' - Subtract gamma from New background
%                       {true|false}. Default is true.
%            'SatPixR'- A map indicating saturated pixels in the Ref.
%                       This can be an image of 0 for unsaturated pixels
%                       and ~0 for saturated pixels.
%                       Default is 0 (scalar) meaning that all the pixels
%                       are unsaturated.
%            'SatPixN'- A map indicating saturated pixels in the New.
%                       This can be an image of 0 for unsaturated pixels
%                       and ~0 for saturated pixels.
%                       Default is 0 (scalar) meaning that all the pixels
%                       are unsaturated.
%            'SatLimitR'- Saturation limit for Ref. If empty then don't
%                       use SetPixR map (if exist). If not empty then this
%                       will override the SetPixR option.
%                       Default is empty. Valu is in ADU units.
%            'SatLimitN'- Saturation limit for New. If empty then don't
%                       use SetPixN map (if exist). If not empty then this
%                       will override the SetPixN option.
%                       Default is empty. Valu is in ADU units.
%            'Eps'    - Small factor to add to the denominator in order
%                       to avoid divsion by zero. Default is 1e-8.
%                       It seems that this factor is not needed, but
%                       it doesn't effect the results.
%            'SigmaPosX' - The astrometric errors in X.
%                       This is used for the calculation of S_corr.
%                       If empty then do not include astrometric errors
%                       in S_corr. Default is empty.
%            'SigmaPosY' - The astrometric errors in Y.
%                       This is used for the calculation of S_corr.
%                       If empty then do not include astrometric errors
%                       in S_corr. Default is empty.
%            'NormS'  - Normalize S to have std of 1. Default is false.
%            'NormScorr' - Normalize Scorr to have std of 1. Default is false.
%            'NormD'  - Normalize D to have std of 1. Default is false.
%            'DelEdge' - Assign some constant to edges of the final
%                       image subtraction products {true|false}.
%                       Default is true.
%            'EdgeBuffer' - Edge buffer width. Default is 15 pixels.
%                       It is recomended to use buffer size of 2-3
%                       PSF widths.
%            'EdgeVal' - Constant value to assign to edges. Default is 0.
%            'Verbose' - Verbose {true|false}. Default is true.
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
%            image2sim.m, images2sim.m, ...
% Output : - A structure containing the following fields:
%            .D       - The proper image subtraction statistics (D).
%            .S       - The transient source detection statistics (S/F_r).
%            .Scorr   - The transient source detection statistics
%                       corrected for source noise and astrometric noise
%                       (S_corr).
%            .fftP_D  - The fft of the PSF of D (P_D).
%            .fftP_Dn - The fft of the PSF of a delta function in N
%                       as appears in D (P_Dn).
%            .fftP_Dr - The fft of the PSF of a delta function in R
%                       as appears in D (P_Dn).
%            .Beta    - The flux-based zero point ratio: F_n/F_r.
%            .Gamma   - Residual background difference: B_n - B_r.
%            .F_S     - Flux scaling needed in order to calculate fluxes.
%            .SigmaFlux - Map of sigma in flux - this is correct only if
%                       the images were provided in units of electrons
%                       or the gain was supplied.
% License: This program can be used only for non-profit educational
%          and scientific use. For any other use please contact:
%          Yeda, Weizmann Institute of Science.
% Knonw issues: 1. The beta-gamma solver (i.e., gain matching) should be
%               improved. Specifically, when the transperancy is very bad
%               (i.e., beta is low) there is a distrurbing correlation
%               between beta and gamma.
%               2. psf_builder.m should be improved.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Sep 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Sub]=imsub_fft('ImBS.fits','ImGS.fits','PSF_N',NewPSF,'PSF_R',RefPSF);
% Reliable: 2
%--------------------------------------------------------------------------


ImageField  = 'Im';              % SIM image field
%HeaderField = 'Header';
%FileField   = 'ImageFileName';
%MaskField   = 'Mask';
BackImField = 'BackIm';
ErrImField  = 'ErrIm';

%--------------------------
%--- Default parameters ---
%--------------------------
DefV.GainN             = 'GAIN';    % New image CCD gain (or header keyword)
DefV.GainR             = 'GAIN';    % Ref image CCD gain (or header keyword)
DefV.Align             = false;     % Align images prior to subtraction
DefV.SubBackN          = true;      % Subtract back from New
DefV.SubBackR          = true;      % Subtract back from Ref
DefV.BackMethod        = 'mode_fit';
DefV.SigmaN            = [];
DefV.SigmaR            = [];
DefV.MethodStD         = 'fit';
DefV.RN_N              = 5;         % Readout noise of New [electrons]
DefV.RN_R              = 5;         % Readout noise of Ref [electrons]
DefV.PSF_N             = [];  
DefV.PSF_R             = [];
DefV.PSFbuildPar       = {'BoxHalfSize',[15 15],'Rad0',15,'ZeroEdgesThresh',2,'ZeroEdgesMethod','rad'};
DefV.BetaMethod        = 'robustfit';   % robustfit | ransac
DefV.Beta              = [];         % gain match to multiply R in order to get N
DefV.BetaTol           = 1e-4;
DefV.RF_Wfun           = 'welsch'; %{}; %'fair';     % robust fit Wfun parameter - set to {} in order to use default.
DefV.MaxIterBeta       = 10;      % maximum number of iterations for beta solver
DefV.UpdateGamma       = true;    % subtrcat updated background from new
DefV.SatPixR           = 0;       % Map of saturated pixels in R
DefV.SatPixN           = 0;       % Map of saturated pixels in N
DefV.SatLimitR         = [];      % Saturation limit of R [ADU]
DefV.SatLimitN         = [];      % Saturation limit of N [ADU]
DefV.Eps               = 1e-8;    % constant for avoiding division by zero - not really necessery
DefV.SigmaPosX         = [];      % X-axis astrometric noise [pix]
DefV.SigmaPosY         = [];      % Y-axis astrometric noise [pix]
DefV.NormS             = false;
DefV.NormScorr         = false;
DefV.NormD             = false;
DefV.DelEdge           = true;    % Set sub images edges to EdgeVal
DefV.EdgeBuffer        = 15;      % Edge buffer [pix]
DefV.EdgeVal           = 0;       % Edge value
DefV.Verbose           = true;
DefV.PreFilter         = []; %'highpass_fermi';  % use [] if no filter - doesn't work yet
DefV.PreFilterPar      = [10 3];
% Read parameters and set default values
% All input arguments are stored in the InPar structure
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

%-------------
% Read images 
%-------------
% Read single reference image
SimRef = image2sim(Ref);
% Read (multiple) new images
SimNew = images2sim(New);
% Number of new images
Nim    = numel(SimNew);

%-------------------------------
% treat saturated pixels in Ref 
%-------------------------------
% (before gain correction) - units are still ADU
if (~isempty(InPar.SatLimitR))
    InPar.SatPixR = double(SimRef.(ImageField)>InPar.SatLimitR);
end


%----------------------
% correct for CCD Gain
%----------------------
% set images gain to 1 - i.e., images are in electron units.
if (~isempty(InPar.GainR))
    [~,SimRef] = sim_gain(SimRef,'Gain',InPar.GainR,'CorrectGain',true,'CorrGainBack',true);
end
if (~isempty(InPar.GainN))
    [GainN,SimNew] = sim_gain(SimNew,'Gain',InPar.GainN,'CorrectGain',true,'CorrGainBack',true);
else
    GainN = 1;
end


%--------------
% align images
%--------------
% sim_align_shift assumes an arbitrary shift and small (< few deg) rotation
if (InPar.Align)
    SimNew = sim_align_shift(SimNew,SimRef,varargin{:});
end

%-----------------
% calculate noise
%-----------------
% noise must be scalar for the entire image
if (isempty(InPar.SigmaN))
    SimNew = sim_std(SimNew,varargin{:},'MethodStD',InPar.MethodStD);
    InPar.SigmaN = [SimNew.(ErrImField)];
end
if (isempty(InPar.SigmaR))
    SimRef = sim_std(SimRef,varargin{:},'MethodStD',InPar.MethodStD);
    InPar.SigmaR = [SimRef.(ErrImField)];
end

% properties of Reference image
% estimate and subtract background
if (isfield_notempty(SimRef,BackImField) && ~InPar.SubBackR),
    % use user supplied background
    SimRefBS.(ImageField) = SimRef.(ImageField) - SimRef.(BackImField);
else
    % calculate background
    SimRefBS = sim_background(SimRef,varargin{:},'SubBack',true,'BackMethod',InPar.BackMethod);
end

% properties of New image
% estimate and subtract background
if (all(isfield_notempty(SimNew,BackImField)) && ~InPar.SubBackN),
    % use user supplied background
    SimNewBS.(ImageField) = SimNew.(ImageField) - SimNew.(BackImField);
else
    % calculate background
    SimNewBS = sim_background(SimNew,varargin{:},'SubBack',true,'BackMethod',InPar.BackMethod);
end

%----------
% get PSFs
%----------
% I am not happy with the performences of psf_builder.m
if (isempty(InPar.PSF_N)),
    PSF = psf_builder(SimNew,InPar.PSFbuildPar{:});
    InPar.PSF_N = {PSF.PSF};
else
    if (isnumeric(InPar.PSF_N)),
        InPar.PSF_N  = {InPar.PSF_N};
    end
end
if (isempty(InPar.PSF_R)),
    PSF = psf_builder(SimRef,InPar.PSFbuildPar{:});
    InPar.PSF_R = {PSF.PSF};
else
    if (isnumeric(InPar.PSF_R)),
        InPar.PSF_R  = {InPar.PSF_R};
    end
end


% there is only one reference image
% rename it and get it size
Ref   = SimRef(1).(ImageField);
RefBS = SimRefBS(1).(ImageField);
SizeR = size(Ref);
% the Ref and New must have the same size
Sy   = SizeR(1);
Sx   = SizeR(2);    
% assume that the PSF size is an odd number
% and pad the New/Ref by 0 or 1 pixels...
PadY = double(Util.array.is_evenint(Sy));
PadX = double(Util.array.is_evenint(Sx));


%---------------------------------
% Subtraction: for each new image
%---------------------------------
% initialize output
Sub = struct('D',cell(Nim,1),...
             'S',cell(Nim,1),...
             'Scorr',cell(Nim,1),...
             'Snoise',cell(Nim,1),...
             'F_S',cell(Nim,1),...
             'SigmaFlux',cell(Nim,1),...
             'fftP_D',cell(Nim,1),...
             'fftP_Dn',cell(Nim,1),...
             'fftP_Dr',cell(Nim,1),...
             'Pn',cell(Nim,1),...
             'Pr',cell(Nim,1),...
             'Beta',cell(Nim,1),...
             'Gamma',cell(Nim,1),...
             'IterBeta',cell(Nim,1),...
             'SigmaN',cell(Nim,1),...
             'SigmaR',cell(Nim,1));
for Iim=1:1:Nim,
    % for each New image
    if (InPar.Verbose),
        fprintf('Subtract Image %d out of %d\n',Iim,Nim);
    end
    IsigN    = min(Iim,numel(InPar.SigmaN)); % The index of InPar.SigmaN
    
    %-------------------------------
    % treat saturated pixels in New
    %-------------------------------
    % units are electrons
    if (~isempty(InPar.SatLimitN)),
        InPar.SatPixN = double(SimNew(Iim).(ImageField)>(InPar.SatLimitN.*GainN(Iim)));
    end    
    
    % rename New and get it size
    New   = SimNew(Iim).(ImageField);
    NewBS = SimNewBS(Iim).(ImageField);
    SizeN = size(New);
    
    %-------------------
    %--- PSF padding ---
    %-------------------
    % THIS SHOULD BE REPLACED WITH THE SHIFT OPERATOR...
    % PSF_N and PSF_R should have the same size
    if all(size(InPar.PSF_N{Iim})==size(InPar.PSF_R{1})),
        % same size
        PSF_N = InPar.PSF_N{Iim};
        PSF_R = InPar.PSF_R{1};
    else
        % not the same size
        error('Can not handle PSF with differnt sizes');
    end
    
    if (all(size(PSF_N)==size(New)) && all(size(PSF_R)==size(Ref))),
        % PSF and images have the same size
        % do nothing
        PadX = 0;
        PadY = 0;
    else
        New   = padarray(New,[PadY PadX],'post');
        NewBS = padarray(NewBS,[PadY PadX],'post');
        SizeN = size(New);
        Ref   = padarray(Ref,[PadY PadX],'post');
        RefBS = padarray(RefBS,[PadY PadX],'post');
        SizeR = size(Ref);

        % shift filter such that it will not introduce a shift to the image
        SizeFilter = size(PSF_R);
        PadImY = 0.5.*(SizeN(1)-(SizeFilter(1)-1)-1);
        PadImX = 0.5.*(SizeN(2)-(SizeFilter(2)-1)-1);
        PSF_R  = padarray(PSF_R,[PadImY PadImX],'both');
        %PSF_R = InPar.FiltNorm.*PSF_R./sum(PSF_R(:));
        PSF_N  = padarray(PSF_N,[PadImY PadImX],'both');
        %PSF_N = InPar.FiltNorm.*PSF_N./sum(PSF_N(:));
        % now Ref, New, PSF_R, PSF_N have the same size
    end
 
    
    % shift the PSF to the edges (origin) of the image
    % so it will not introduce a shift
    PSF_N   = fftshift(fftshift(PSF_N,1),2);  
    PSF_R   = fftshift(fftshift(PSF_R,1),2);
    
    % renormalize PSF (just in case)
    PSF_N    = PSF_N./sum(PSF_N(:));
    PSF_R    = PSF_R./sum(PSF_R(:));

    % fft2 of components
    FPn   = fft2(PSF_N);    % FFT of P_n
    FPr   = fft2(PSF_R);    % FFT of P_r
    FN    = fft2(NewBS);    % FFT of N (back. subtracted)
    FR    = fft2(RefBS);    % FFT of R (back. subtracted)
    FNnbs = fft2(New);      % FFT of N (not back. sub.)
    FRnbs = fft2(Ref);      % FFT of R (not back. sub.)
    
    % Pre filtering
    % analog for background subtraction
    % doesn't work yet - low priority
    if (~isempty(InPar.PreFilter)),
        Filter=filter2_design(SizeN,InPar.PreFilter,InPar.PreFilterPar);
        warning('BUG - REQUIRE fftshift!!!!!!');
        FN = FN.*Filter;
        FR = FR.*Filter;
    end
    
    FPrFN = FPr.*FN;   % fft(P_r).*fft(N)
    FPnFR = FPn.*FR;   % fft(P_n).*fft(R)

    %-----------------------------------
    %--- gain matching (zero points) ---
    %-----------------------------------
    % find beta and gamma
    if (isempty(InPar.Beta)),
        % assume gain matching Beta factor is unknown and try to fit it
        switch lower(InPar.BetaMethod)
            case {'robustfit','ransac'}
                % solve for beta and gamma (no translations)
                % using iterative robust fitting
                
                % set initial guess for the value of beta and gamma
                Beta     = 1;
                GammaTag = 0;
                
                BetaGammaConverg = false;
                IterBeta = 0;
                while (~BetaGammaConverg && IterBeta<=InPar.MaxIterBeta),
                    % The denominator of D with the beta term:
                    IterBeta = IterBeta + 1;
                    Den     = Beta.^2.*InPar.SigmaR.^2 .* FPn.*conj(FPn) + InPar.SigmaN(IsigN).^2 .* FPr.*conj(FPr) + InPar.Eps;
                    SqrtDen = sqrt(Den);
                    Dn = ifft2(FPrFN./SqrtDen);
                    Dr = ifft2(FPnFR./SqrtDen);

                    % select only pixels which are >1 sigma above the image
                    % noise.
                    % Is this necesssery?
                    [Back_Dn,Std_Dn] = mode_fit(Dn);
                    [Back_Dr,Std_Dr] = mode_fit(Dr);
                    %Il0 = Dn(:)>1 & Dr(:)>1 & InPar.SatPixR(:)==0 & InPar.SatPixN(:)==0;
                    Il0 = Dn(:)>Std_Dn.*7 & Dr(:)>Std_Dr.*7 & InPar.SatPixR(:)==0 & InPar.SatPixN(:)==0;
                    %Il0 = Dn(:)>InPar.SigmaN(IsigN) & Dr(:)>InPar.SigmaR & InPar.SatPixR(:)==0 & InPar.SatPixN(:)==0;
                    %Il0 = Dn(:)>-Inf & Dr(:)>-Inf;
                    
                    % robust fitting of the Dn and Dr
                    
                    switch lower(InPar.BetaMethod)
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
                    
%                      plot(Dr(Il0),Dn(Il0),'.');
%                      hold on
%                      plot(Dr(Il0),Par(1) + Par(2).*Dr(Il0),'r.')
%                      input('next')
                end
                Gamma = GammaTag.*sqrt(InPar.SigmaN(IsigN).^2 + Beta.^2.*InPar.SigmaR.^2);
                % update the background level of the New
                if (InPar.UpdateGamma),
                    FN    = fft2(NewBS - Gamma);
                    FPrFN = FPr.*FN;
                end
                
                %I10 = find(abs(Par(1) + Par(2).*Dr(Il0)-Dn(Il0))<10);
                %Par = robustfit(Dr(Il0(I10)),Dn(Il0(I10)),InPar.RF_Wfun);
                %GammaTag = Par(1);
                %Beta     = Par(2);
                    
                
                plot(Dr(Il0),Dn(Il0),'.');
                hold on
                plot(Dr(Il0),Par(1) + Par(2).*Dr(Il0),'r.')
                
                
                if (InPar.Verbose),
                    fprintf('Gain matching using %s method\n',InPar.BetaMethod);
                    fprintf('   Solution found after %d iterations\n',IterBeta);
                    fprintf('   gamma_tag (offset) = %f\n',Par(1));
                    fprintf('   beta (N/R gain)    = %f\n',Par(2));
                end
                
            otherwise
                error('Unknown BetaMethod option');
        end
    else
        Beta = InPar.Beta;
        Gamma = 0;
        IterBeta = 0;
        
    end
    
    
    %---------------------------------------------------------------------
    %--- Calculate the Score image statistics (S) for source detection ---
    %---------------------------------------------------------------------
    FPn2 = FPn.*conj(FPn);       % |fft(P_n)|^2
    FPr2 = FPr.*conj(FPr);       % |fft(P_r)|^2
    Den     = Beta.^2.*InPar.SigmaR.^2 .* FPn2 + InPar.SigmaN(IsigN).^2 .* FPr2 + InPar.Eps;
    SqrtDen = sqrt(Den);         % don't use 1/sqrt() as this will introduce division by zero...
    Kn  = conj(FPn).*FPr2./Den;  % This is not the Kn definition in the paper (missing Fn,Fr)
    Kr  = conj(FPr).*FPn2./Den;  % This is not the Kr definition in the paper (missing Fn,Fr)

    % The filtered image subtraction statistics S:
    % This is the log-likelihood to find a transient (with the image PSF)
    % at each pixel.
    % alternatively can calculated by filtering D with P_d:
    % Sub(Iim).S = ifft2(Kn.*FN - Beta.*Kr.*FR);   % formally this is S/F_n
    Sub(Iim).S = ifft2(Beta.*Kn.*FN - Beta.^2.*Kr.*FR);   % formally this is S/F_r
    
    
    %------------------------------------------------------
    %--- correct for source noise and astrometric noise ---
    %------------------------------------------------------
    % This assumes that at this stage the images are in units of electrons
    F_r = 1;
    F_n = F_r.*Beta;
    Kr2 = fft2(ifft2(F_r.*Beta.^2.*Kr).^2);   % <----------
    Kn2 = fft2(ifft2(F_n         .*Kn).^2);   % <----------
    %Kr2 = fft2(ifft2(Kr).^2);
    %Kn2 = fft2(ifft2(Kn).^2);
    
    
    %Kn2 = Kn2./sum(Kn2(:));
    %Kr2 = Kr2./sum(Kr2(:));
    
    
    if (isempty(InPar.SigmaPosX) || isempty(InPar.SigmaPosY)),
        %Vadd = 0;
        FgradN = 0;
        FgradR = 0;
    else
        IsigposX = min(Iim, numel(InPar.SigmaPosX));
        IsigposY = min(Iim, numel(InPar.SigmaPosY));
        
        % Kn is in Fourier space - the gradient is calculated in real space
        [GradNx,GradNy] = gradient(ifft2(F_n         .*Kn.*FN));  % gradient of S_n
        [GradRx,GradRy] = gradient(ifft2(F_r.*Beta.^2.*Kr.*FR));  % gradient of S_r
        %[GradNx,GradNy] = gradient(ifft2(Kn.*FN));  % gradient of S_n      <-----------------
        %[GradRx,GradRy] = gradient(ifft2(Kr.*FR));  % gradient of S_r      <-----------------
        FgradN = (InPar.SigmaPosX(IsigposX).*GradNx).^2 + (InPar.SigmaPosY(IsigposY).*GradNy).^2;
        FgradR = (InPar.SigmaPosX(IsigposX).*GradRx).^2 + (InPar.SigmaPosY(IsigposY).*GradRy).^2;
       
        %VaddFT = Kn2.*FgradN + Kr2.*FgradR;
        % Any additional noise (e.g., color-refraction noise) should be
        % added to Vadd.
        %Vadd = FgradN + FgradR;
        
    end
    
    % Variance correction
    % Incluing: Source noise and astrometric noise
    %Sub(Iim).Snoise = sqrt(ifft2(Kn2.*FNnbs + Kr2.*FRnbs) + Vadd);
    Sub(Iim).Snoise = sqrt(ifft2(Kn2.*FNnbs + Kr2.*FRnbs) + FgradN + FgradR);
    
    % estimate the ratio between N and R noise:
    %sqrt(mediannd((ifft2(Kn2.*FNnbs) + FgradN)./(ifft2(Kr2.*FRnbs) + FgradR)))
    
    
    %---------------------------
    % Transient flux estimation
    %---------------------------
    % Flux normalization (F_S):
    % the abs value is needed in order to get rid of small imaginary part
    Sub(Iim).F_S = abs(sumnd((Beta.^2.* FPn2 .* FPr2)./Den));
    
    % Flux error
    % assumes images are in units of electrons or gain is known
    Sub(Iim).SigmaFlux = Sub(Iim).Snoise./Sub(Iim).F_S;
    
    
    %----------------------------------------
    %--- Image subtraction statistics (D) ---
    %----------------------------------------
    % The proper image subtraction statistics which have
    % uncorelated noise in the background-dominated noise limit.
    Sub(Iim).D      = ifft2((FPrFN - Beta.*FPnFR)./SqrtDen);
    
    %----------------------
    %--- remove padding ---
    %----------------------
    Sub(Iim).S      = Sub(Iim).S(1:end-PadY,1:1:end-PadX);    
    Sub(Iim).Snoise = Sub(Iim).Snoise(1:end-PadY,1:1:end-PadX);
    Sub(Iim).D      = Sub(Iim).D(1:end-PadY,1:1:end-PadX);
    Sub(Iim).Scorr  = Sub(Iim).S./Sub(Iim).Snoise;
    
    warning('PSF padding is not corrected yet!')
    % specifically, the PSF may be larger by one pixel compared with
    % the image.
    % This will be fixed with the new PSF class.
    
    % calculate the PSFs of ImSub
    % PSF for unresolved sources
    Sub(Iim).fftP_D   = FPn.*FPr./SqrtDen;
    % PSF for delta function in N
    Sub(Iim).fftP_Dn  = FPr./SqrtDen;
    % PSF for delta function in R
    Sub(Iim).fftP_Dr  = FPn./SqrtDen;
    
    % normalize output std
    if (InPar.NormS),
         [Back,Std] = Util.stat.mode_fit(Sub(Iim).S);
         Sub(Iim).S = (Sub(Iim).S - Back)./Std;
    end
    if (InPar.NormScorr),
         [Back,Std] = Util.stat.mode_fit(Sub(Iim).Scorr);
         Sub(Iim).Scorr = (Sub(Iim).Scorr - Back)./Std;
    end
    if (InPar.NormD),
         [Back,Std] = Util.stat.mode_fit(Sub(Iim).D);
         Sub(Iim).D = (Sub(Iim).D - Back)./Std;
    end
    
    % set edges to constant (to avoid convolution edge effects)
    if (InPar.DelEdge),
        SizeD = size(Sub(Iim).D);
        [MatX,MatY] = meshgrid((1:1:SizeD(2)),(1:1:SizeD(1)));
        FlagEdge = MatX>InPar.EdgeBuffer & MatX<(SizeD(2)-InPar.EdgeBuffer) & ...
                   MatY>InPar.EdgeBuffer & MatY<(SizeD(1)-InPar.EdgeBuffer);
        Sub(Iim).D(~FlagEdge)     = InPar.EdgeVal;
        Sub(Iim).S(~FlagEdge)     = InPar.EdgeVal;
        Sub(Iim).Scorr(~FlagEdge) = InPar.EdgeVal;
    end
    
    % Store additional info:
    Sub(Iim).Beta     = Beta;
    Sub(Iim).Gamma    = Gamma;
    Sub(Iim).Pn       = fftshift(PSF_N);
    Sub(Iim).Pr       = fftshift(PSF_R);
    Sub(Iim).IterBeta = IterBeta;
    Sub(Iim).SigmaR   = InPar.SigmaR;
    Sub(Iim).SigmaN   = InPar.SigmaN(IsigN);
    
    % warning flgas
%     if (abs(Sub(Iim).Gamma)<(3.*Sub(Iim).SigmaN)),
%         % Gamma is very different than initial value
%     end
    
end
    
    
    

