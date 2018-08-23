function Res=fit_sn_rise(LC,varargin)
%--------------------------------------------------------------------------
% fit_sn_rise function                                              FitFun
% Description: Fit various functions appropriate for the rise of SN
%              light curve. The fitted functions are:
%              L_max*(1-((t-t_max)/t_rise)^2)
%              L_max.*(1-exp(-(t-t_start)./t_rise))
%              L_max.*erfc(sqrt(t_diff./(2.*(t-t_start))))
% Input  : - Supernova light curve [time(day), mag, err]
%          * Arbitrary number of pair of arguments, ...,key,val,...
%            The following keywords are available:
%            'Type'  - LC type {'mag','lum'}. Default is 'mag'.
%                      If 'lum', then the first inpt arguments is
%                      [time(day), Lum(erg/s), err(erg/s)]
%            'DM'    - Distance modulus [mag] to use in conversion from
%                      magnitude to luminosity. Default is 0.
%            'A'     - Extinction in band [mag] to use in conversion from
%                      magnitude to luminosity. Default is 0.
%            'SunAbsMag' - Sun abs. mag. in band. Default is 4.66 (PTF R).
%            'Nsim'  - Number of bootstrap simulations for error
%                      estimation. Default is 100.
%            'Options' - optimset structure to pass to the non linear
%                        fitting.
% Output : - A structure with the following fields:
%            .t2   - a t^2 law fit
%            .exp  - an exp fit
%            .erfc - an erfc fit
%            Each field contains the following seb fields:
%            .Par     - Best fit parameters
%                       in .t2 these are [L_max, t_max, t_rise]
%                       in .exp these are [L_max, t_start, t_rise]
%                       in .erfc these are [L_max, t_diff, t_start]
%            .ParErrB - Bootstrap errors in best fit parameters
%            .Chi2    - \chi^2 of best fit
%            .Dof     - number of degrees of freedom
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Jul 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=fit_sn_rise(LC);
% Reliable: 2
%--------------------------------------------------------------------------

DefV.Type      = 'mag';
DefV.DM        = 0;
DefV.A         = 0;
DefV.SunAbsMag = 4.66;
DefV.Nsim      = 100;
DefV.Options   = [];
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

switch lower(InPar.Type)
    case 'mag'     
        Time   = LC(:,1);
        Mag    = LC(:,2);
        MagErr = LC(:,3);

        SolarLum = 3.84e33;

        Flux    = 10.^(-0.4.*(Mag - InPar.DM - InPar.A - InPar.SunAbsMag)).*SolarLum;
        FluxErr = MagErr.*Flux;
    case 'lum'
        Time    = LC(:,1);
        Flux    = LC(:,2);
        FluxErr = LC(:,3);
    otherwise
        error('Unknown Type option');
end
N       = length(Time);

       
%--- fit a L_max*(1-((t-t_max)/t_rise)^2) law ---
% Par(1) = L_max
% Par(2) = t_max
% Par(3) = t_rise
Fun     = @(Par,t) Par(1).*(1-((t-Par(2))./Par(3)).^2);
FunChi2 = @(Par,t,Flux,FluxErr) sum(((Flux-Par(1).*(1-((t-Par(2))./Par(3)).^2))./FluxErr).^2);
[MaxLum, MaxLumInd] = max(Flux);
GuessPar(1) = MaxLum;                    % guess L_max
GuessPar(2) = Time(MaxLumInd);           % guess t_max
GuessPar(3) = GuessPar(2)-min(Time) + 5; % guess t_rise

Data = [Time, Flux, FluxErr];

[Beta,Chi2,E,O]=fminsearch_chi2(Time,Flux,FluxErr,{Fun},GuessPar,InPar.Options);

% bootstrap errors
Theta = zeros(InPar.Nsim,length(Beta));
for Isim=1:1:InPar.Nsim,
     Rand    = rand(N,1);
     NewInd  = ceil(Rand.*N);
     NewData = Data(NewInd,:);
    [BetaS,Chi2S,ES,OS]=fminsearch_chi2(NewData(:,1),NewData(:,2),NewData(:,3),{Fun},GuessPar,InPar.Options);

    Theta(Isim,:) = BetaS;
end
MeanTheta = mean(Theta);
StD       = sqrt(sum((bsxfun(@minus,Theta,MeanTheta)).^2)./(InPar.Nsim-1));

StDrob(1)    = diff(err_cl(Theta(:,1)-MeanTheta(1),0.6827),1,2)./2;
StDrob(2)    = diff(err_cl(Theta(:,2)-MeanTheta(2),0.6827),1,2)./2;
StDrob(3)    = diff(err_cl(Theta(:,3)-MeanTheta(3),0.6827),1,2)./2;

t0    = Beta(2) - Beta(3);
t0Err = sqrt(StD(2).^2 + StD(3).^2);
% calculate the time it took the LC to rise by factor of 1/e relative to
% max light
Alpha = 1./exp(1);
t_e    = Beta(2) - Beta(3).*sqrt(1-Alpha) - t0;
t_eErr = sqrt( StD(2).^2 + (StD(3).*sqrt(1-Alpha)).^2);

% store results
Res.t2.Par     = Beta;   % [L_max, t_max, t_rise] 
Res.t2.ParErrB = StD;
Res.t2.ParErrBr= StDrob;
Res.t2.Chi2    = Chi2;
Res.t2.Dof     = N - length(Beta);
Res.t2.Fun     = Fun;
Res.t2.t0      = t0;
Res.t2.t0Err   = t0Err;
Res.t2.t_e     = t_e;
Res.t2.t_eErr  = t_eErr;



%--- fit exponential ---
% F = L_max.*(1-exp(-(t-t_start)./t_rise));
% Par(1) = L_max
% Par(2) = t_start
% Par(3) = t_rise

Fun     = @(Par,t) Par(1).*(1-exp(-(t-Par(2))./Par(3)));
FunChi2 = @(Par,t,Flux,FluxErr) sum(((Flux- Par(1).*(1-exp(-(t-Par(2))./Par(3))) )./FluxErr).^2);
[MaxLum, MaxLumInd] = max(Flux);
[MinLum, MinLumInd] = min(Flux);
GuessPar(1) = MaxLum;                      % guess L_max
GuessPar(2) = Time(MinLumInd)- 5;         % guess t_start
GuessPar(3) = (Time(MaxLumInd)-GuessPar(2))./3; % guess t_rise


Data = [Time, Flux, FluxErr];

[Beta,Chi2,E,O]=fminsearch_chi2(Time,Flux,FluxErr,{Fun},GuessPar,InPar.Options);

% bootstrap errors
Theta = zeros(InPar.Nsim,length(Beta));
for Isim=1:1:InPar.Nsim,
     Rand    = rand(N,1);
     NewInd  = ceil(Rand.*N);
     NewData = Data(NewInd,:);
    [BetaS,Chi2S,ES,OS]=fminsearch_chi2(NewData(:,1),NewData(:,2),NewData(:,3),{Fun},GuessPar,InPar.Options);

    Theta(Isim,:) = BetaS;
end
MeanTheta = mean(Theta);
StD       = sqrt(sum((bsxfun(@minus,Theta,MeanTheta)).^2)./(InPar.Nsim-1));

StDrob(1)    = diff(err_cl(Theta(:,1)-MeanTheta(1),0.6827),1,2)./2;
StDrob(2)    = diff(err_cl(Theta(:,2)-MeanTheta(2),0.6827),1,2)./2;
StDrob(3)    = diff(err_cl(Theta(:,3)-MeanTheta(3),0.6827),1,2)./2;

t0        = Beta(2);
t0Err     = StD(2);

% store results
Res.exp.Par     = Beta;   % [L_max, t_max, t_rise] 
Res.exp.ParErrB = StD;
Res.exp.ParErrBr= StDrob;
Res.exp.Chi2    = Chi2;
Res.exp.Dof     = N - length(Beta);
Res.exp.Fun     = Fun;
Res.exp.t0      = t0;
Res.exp.t0Err   = t0Err;

%--- fit L_max.*erfc(sqrt(t_diff./(2.*(t-t_start)))) ---
% F = L_max.*erfc(abs(sqrt(t_diff./(2.*(t-t_start))))); % Piro & Nakar (2013)
% Par(1) = L_max
% Par(2) = t_diff
% Par(3) = t_start

Fun     = @(Par,t) Par(1).*erfc(abs(sqrt(Par(2)./(2.*(t-Par(3))))));
FunChi2 = @(Par,t,Flux,FluxErr) sum(((Flux-  Par(1).*erfc(abs(sqrt(Par(2)./(2.*(t-Par(3))))))  )./FluxErr).^2);
[MaxLum, MaxLumInd] = max(Flux);
[MinLum, MinLumInd] = min(Flux);
GuessPar(1) = MaxLum;                      % guess L_max
GuessPar(2) = (max(Time) - min(Time))./10; % guess t_diff
GuessPar(3) = min(Time) - 5;               % guess t_start

Data = [Time, Flux, FluxErr];

[Beta,Chi2,E,O]=fminsearch_chi2(Time,Flux,FluxErr,{Fun},GuessPar,InPar.Options);

% bootstrap errors
Theta = zeros(InPar.Nsim,length(Beta));
for Isim=1:1:InPar.Nsim,
     Rand    = rand(N,1);
     NewInd  = ceil(Rand.*N);
     NewData = Data(NewInd,:);
    [BetaS,Chi2S,ES,OS]=fminsearch_chi2(NewData(:,1),NewData(:,2),NewData(:,3),{Fun},GuessPar,InPar.Options);

    Theta(Isim,:) = BetaS;
end
MeanTheta = mean(Theta);
StD       = sqrt(sum((bsxfun(@minus,Theta,MeanTheta)).^2)./(InPar.Nsim-1));


t0        = Beta(3);
t0Err     = StD(3);

% store results
Res.erfc.Par     = Beta;   % [L_max, t_max, t_rise] 
Res.erfc.ParErrB = StD;
Res.erfc.Chi2    = Chi2;
Res.erfc.Dof     = N - length(Beta);
Res.erfc.Fun     = Fun;
Res.erfc.t0      = t0;
Res.erfc.t0Err   = t0Err;







