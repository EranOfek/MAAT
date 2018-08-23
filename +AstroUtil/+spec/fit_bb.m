function [Res,SpecBB]=fit_bb(Spec,varargin)
% Fit a blackbody to spectrum
% Package: AstroUtil.spec
% Description: Fit a black body spectrum to a list of spectral measurments,
%              spectrum or photometric measurments.
% Input  : - The spectral points to fit [Wavelength, Flux, [Err]].
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'FitMethod' - One of the following fit methods:
%                          'lscov' - use lscov. Default.
%                          '\'
%                          'ratio_chi2'
%                          'ratio_rms'
%                          'ratio_relrms'
%            'BadRanges' - Two column matrix of bad ranges to remove from
%                          spectrum before fitting. Line per bad range.
%                          This can be used to mask lines.
%                          Default is [].
%            'WaveUnits' - Wavelength units. Default is 'Ang'.
%                          See convert.units.m for options.
%            'IntUnits'  - Flux units. Default is 'erg*cm^-2*s^-1*Ang^-1'.
%                          See convert.flux.m for options.
%                          E.g., 'AB', 'STmag','mJy','ph/A',...
%            'Trange'    - Temperature range [K] in which to search for
%                          solution. Default is [100 5e6].
%            'Tpoint'    - Number of points in each search range.
%                          Default is 5.
%            'Thresh'    - Convergence threshold. Default is 0.01;
%            'ColW'      - Column index in spectrum input containing
%                          the wavelength. Default is 1.
%            'ColF'      - Column index in spectrum input containing
%                          the flux. Default is 2.
%            'ColE'      - Column index in spectrum input containing
%                          the error in flux. Default is 3.
%            'RelErrFlux' - If the error column is empty then will se the
%                          error to be the flux multiplied by this factor.
%                          If this factor is empty then will set the errors
%                          to 1. Default is 0.05.
%            'FitFun'    - Function for ratio calculation in initial
%                          iterative search. Default is @mean.
%            'Nsigma'    - Return the errors for Nsigma confidence
%                          interval. Default is 1.
%            'Plot'      - plot specectra and best fit. Default is false.
% Output : - A structure containing the best fit black body results.
%          - AstSpec class object containing the calculated best fit
%            black-body spectrum.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res,SpecBB]=AstroUtil.spec.fit_bb([AS(1).Wave,AS(1).Int]);
%          plot(AS(1)); hold on; plot(SpecBB);
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.FitMethod          = 'lscov';
DefV.BadRanges          = [];
DefV.WaveUnits          = 'Ang'; 
DefV.IntUnits           = 'erg*cm^-2*s^-1*Ang^-1';
DefV.Trange             = [100 5e6];
DefV.Tpoint             = 5;
DefV.Thresh             = 0.01;
DefV.ColW               = 1;
DefV.ColF               = 2;
DefV.ColE               = 3;
DefV.RelErrFlux         = []; %0.05;   % relative error of flux to use if error not provided or empty
DefV.FitFun             = @mean;  % @median
DefV.Nsigma             = 1;      % number of sigmas for error calculation
DefV.Plot               = false;

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% check if Spec contains Errors - if not add equal weights
Nw  = size(Spec,1);
if (size(Spec,2)==2)
    if (isempty(InPar.RelErrFlux))
        Spec = [Spec, ones(Nw,1)];
    else
        Spec = [Spec, Spec(:,InPar.ColF).*InPar.RelErrFlux];
        Ferr0 = Spec(:,InPar.ColE)==0;
        Spec(Ferr0,InPar.ColE) = max(Spec(Ferr0,InPar.ColE), max(Spec(:,InPar.ColF)).*InPar.RelErrFlux  );
    end
end

% remove bad ranges
IndBad = Util.array.find_ranges(Spec(:,InPar.ColW),InPar.BadRanges);
Spec(IndBad,2) = NaN;
Spec = Spec(~isnan(Spec(:,2)),:);
Nw   = size(Spec,1);
Npar = 2;
Ndof = Nw - Npar;



Cont = true;
VecChi2 = zeros(InPar.Tpoint,1);
Trange  = logspace(log10(InPar.Trange(1).*0.5),log10(InPar.Trange(2).*2),InPar.Tpoint).';
Iter    = 0;
while Cont
    % calculate chi^2 as a function of T in some range
    % this is required for robustness
    Iter = Iter + 1;
    for It=1:1:InPar.Tpoint
        [VecChi2(It)] = do_bb_fit(Spec,Trange(It),InPar);
    end
    
    %[Trange, VecChi2]
    TrangeRatio = Trange(2)./Trange(1);
    [~,MinInd] = min(VecChi2);
    BestT      = Trange(MinInd);
    TrangeOld  = Trange;
    Trange     = logspace(log10(BestT./TrangeRatio),log10(BestT.*TrangeRatio),InPar.Tpoint).';
    
    Cont = (TrangeRatio-1)>InPar.Thresh;
end

% re fit with normzlized errors
switch lower(InPar.FitMethod)
    case {'ratio_chi2','lscov'}
        Chi2 = min(VecChi2);
        Spec(:,InPar.ColE) = Spec(:,InPar.ColE).*sqrt(Chi2./Ndof);
        for It=1:1:InPar.Tpoint
            [VecChi2(It)] = do_bb_fit(Spec,Trange(It),InPar);
        end
        
    otherwise
        % do nothing
end

MeanT = mean(Trange);
Trange = Trange - MeanT;
Par=[Trange.^2, Trange, ones(size(Trange))]\VecChi2;

% fit parabola to find best fit and errors
% MeanT = mean(Trange);
% Par = polyfit(Trange-MeanT,VecChi2,2);
BestT = MeanT-Par(2)./(2.*Par(1));

% final fit at BestT
[BestChi2,Ratio,RatioErr,Resid,BB]=do_bb_fit(Spec,BestT,InPar);

Res.BestT       = BestT;
Res.Chi2        = BestChi2;
Res.Ndof        = Ndof;
Res.Npar        = Npar;
Res.Ratio       = Ratio;
Res.RatioErr    = RatioErr;
Res.AngRad      = sqrt(Ratio).*RAD.*3600;   % [arcsec]
Res.Resid       = Resid;

switch lower(InPar.FitMethod)
    case {'ratio_chi2','lscov'}
        % fit parabola to estimate errors
        ProbSigma = 1-(1-normcdf(InPar.Nsigma,0,1)).*2;
        %Par = polyfit(Trange,VecChi2,2);
        Par1 = Par;
        Par1(3) = Par1(3) - Ndof - chi2cdf(ProbSigma,Npar);
        Troots  = MeanT+roots(Par1);

        Terr  = [BestT-min(Troots), max(Troots)-BestT];

    otherwise
        % can't estimate errors in non-chi2 mode
        Terr = [NaN NaN];
end
Res.Terr = Terr;

if (nargout>1)
    SpecBB               = AstSpec;
    SpecBB.Wave          = Spec(:,InPar.ColW);
    SpecBB.Int           = Ratio.*BB;
    SpecBB.WaveUnits     = InPar.WaveUnits;
    SpecBB.IntUnits      = InPar.IntUnits;
    SpecBB.source        = 'AstSpec/blackbody';
    SpecBB.z             = 0;
end
    
if (InPar.Plot)
    plot(Spec(:,1),Spec(:,2))
    hold on
    plot(Spec(:,1),BB.*Ratio)
end



    
end

    
function [VecChi2,Ratio,RatioErr,Resid,BB]=do_bb_fit(Spec,Temp,InPar)

    BB = AstSpec.blackbody(Temp,Spec(:,InPar.ColW),InPar.IntUnits,InPar.WaveUnits,'mat');
    BB = BB(:,2);
    
    RatioErr = NaN;
    switch lower(InPar.FitMethod)
        case 'ratio_chi2'
            Ratio = InPar.FitFun(Spec(:,InPar.ColF)./BB);
            Resid = Spec(:,InPar.ColF) - Ratio.*BB;
            VecChi2  = sum((Resid./Spec(:,InPar.ColE)).^2);
        case 'ratio_rms'
            Ratio = InPar.FitFun(Spec(:,InPar.ColF)./BB);
            Resid = Spec(:,InPar.ColF) - Ratio.*BB;
            VecChi2  = nanstd(Resid);
        case 'ratio_relrms'
            Ratio = InPar.FitFun(Spec(:,InPar.ColF)./BB);
            Resid = Spec(:,InPar.ColF) - Ratio.*BB;
            Fn0   = Spec(:,InPar.ColF)>0;
            VecChi2  = nanstd(Resid(Fn0)./Spec(Fn0,InPar.ColF));
        case '\'
            Ratio = BB\Spec(:,InPar.ColF);
            Resid = Spec(:,InPar.ColF) - Ratio.*BB;
            VecChi2  = nanstd(Resid);
        case 'lscov'
            [Ratio,RatioErr] = lscov(BB,Spec(:,InPar.ColF),1./(Spec(:,InPar.ColE).^2));
            Resid = Spec(:,InPar.ColF) - Ratio.*BB;
            VecChi2  = sum((Resid./Spec(:,InPar.ColE)).^2);
        otherwise
            error('Unknown FitMethod option');
    end
end