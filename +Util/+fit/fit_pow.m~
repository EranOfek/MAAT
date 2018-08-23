function [Res]=fit_pow(X,Y,Err,varargin)
%--------------------------------------------------------------------------
% fit_pow function                                                  FitFun
% Description: Fit power-law function of the form Y=A*X^Alpha to data
%              using non linear lesat squares.
% Input  : - Vector of X data.
%          - Vector of Y data.
%          - Vector of Err on Y, if scalar, then use the scalar error
%            for all Y values. Default is 1.
%          * Arbitrary number of pairs of optional input arguments,
%            ...,key,val,... The following keywords are available:
%            'Par0'     - Guess parameters [A, Alpha]. Default is [1 1].
%            'ResidType'- Type of residuals to use in the sigma clipping:
%                         'Abs'  - absolute residuals.
%                         'Rel'  - relative residuals (default).
%            'MaxIter'  - Maximum number of sigma clipping iterations.
%                         If 0, then no sigma clipping. Default is 0.
%            'Method'   - Sigma clipping method. See clip_resid.m for
%                         options.
%            'Mean'     - Sigma clipping mean of sample estimator.
%                         See clip_resid.m for options.
%            'Clip'     - Lower and upper clipping values.
%                         See clip_resid.m for options.
%            'StdZ'     - Add epsilon for sigma clipping.
%                         See clip_resid.m for options.
%            'CalcBS'   - Calculate bootstrap errors.
%                         {'y' | 'n'}. Default is 'y'.
%            'Nsim'     - Number of bootstarp simulations. Default is 1000.
% Output : - Structure containing the following fields:
%            .Par     - Best fit parameters [A, Alpha].
%            .ParErr  - Error in best fit parameters
%            .Resid   - Vector of absolute residuals.
%            .RelResid- Vector of relative residuals.
%            .H       - Hessian matrix.
%            .Cov     - Covariance matrix.
%            .Chi2    - \chi^{2}
%            .ExitFlag- extit flag from last fminsearch call.
%            .Output  - Output from last fminsearch call.
%            .Nobs    - Number of observations.
%            .Npar    - Number of free parameters.
%            .Dof     - Number of degrees of freedom.
%            .ParErrBS- Errors in parameters estimated using Bootstrap.
%            .ParBiasBS- Bias in parameters estimated using Bootstrap.
% Tested : Matlab 2011b
%     By : Eran O. Ofek                     May 2013
%    URL : http://weizamann.ac.il/home/eofek/matlab/
% Example: X=logspace(0,2,100)'; Y=5.*X.^0.85+randn(100,1).*0.5;
%          Res=fit_pow(X,Y,0.5);
% Reliable: 2
%--------------------------------------------------------------------------

if (numel(Err)==1),
    Err = Err.*ones(size(X));
end


DefV.Par0      = [1 1];
DefV.ResidType = 'Rel';
DefV.MaxIter   = 0;
DefV.Method    = 'StD';
DefV.Mean      = 'Median';
DefV.Clip      = [Inf Inf];
DefV.StdZ      = 'y';
DefV.CalcBS    = 'y';
DefV.Nsim      = 1000;
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});


Fun = @(Par,X)Par(1).*X.^Par(2);
Chi2Fun = @(Par,X,Y,Err) sum(((Y-Par(1).*X.^Par(2))./Err).^2);
[Beta,Chi2,ExitFlag,Output]=fminsearch_chi2(X,Y,Err,Fun,InPar.Par0);

Res.Par      = Beta;   

Res.Resid  = Y - Res.Par(1).*X.^Res.Par(2);
Res.RelResid = Res.Resid./(Res.Par(1).*X.^Res.Par(2));
switch lower(InPar.ResidType)
    case 'rel'
        Resid = Res.RelResid;
    case 'abs'
        Resid = Res.Resid;
    otherwise
        error('Unknwon ResidType option');
end

% sigma clipping
for Iter=1:1:InPar.MaxIter,
   [IndGR,GR,ResidFlag] = clip_resid(Resid,varargin{:});
   [Beta,Chi2,ExitFlag,Output]=fminsearch_chi2(X(IndGR),Y(IndGR),Err(IndGR),Fun,InPar.Par0);

   Res.Par      = Beta;   
   Res.Resid  = Y - Res.Par(1).*X.^Res.Par(2);
   Res.RelResid = Res.Resid./(Res.Par(1).*X.^Res.Par(2));
   Res.ResidFlag = ResidFlag;
   switch lower(InPar.ResidType)
       case 'rel'
           Resid = Res.RelResid;
       case 'abs'
           Resid = Res.Resid;
       otherwise
           error('Unknwon ResidType option');
   end    
end

% calculate errors
[H]=calc_hessian(Chi2Fun,Beta,[],X,Y,Err);
Res.H = H;
Res.Cov  = inv(0.5.*H);
Res.ParErr = sqrt(diag(Res.Cov));

% additional info
Res.Chi2     = Chi2;
Res.ExitFlag = ExitFlag;
Res.Output   = Output;
     
Res.Nobs     = length(X);
Res.Npar     = length(Res.Par);
Res.Dof      = Res.Nobs - Res.Npar;

% Bootstrap errors
switch lower(InPar.CalcBS)
    case 'y'
       for Isim=1:1:InPar.Nsim,

           Rand    = rand(Res.Nobs,1);
           NewInd  = ceil(Rand.*Res.Nobs);
           
           ResBS(Isim) = fit_pow(X(NewInd),Y(NewInd),Err(NewInd),varargin{:},'CalcBS','n');
           if (size(ResBS(Isim).Par,1)==1),
              ResBS(Isim).Par = ResBS(Isim).Par';
           end
       end
       
       Theta = [ResBS.Par];
       MeanTheta = mean(Theta,2);
       for Ipar=1:1:Res.Npar,
          Res.ParErrBS(Ipar) = sqrt(sum((Theta(Ipar,:) - MeanTheta(Ipar)).^2)./(InPar.Nsim-1));
          Res.ParBiasBS(Ipar) = MeanTheta(Ipar) - Res.Par(Ipar);
       end
       
    otherwise
        % do nothing
end
