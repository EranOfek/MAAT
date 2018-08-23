function [MeanJD,Par,ParErr,Chi2,Dof,Resid]=fit_pm(JD,RA,Dec,ErrRA,ErrDec);
%------------------------------------------------------------------------------
% fit_pm function                                                       FitFun
% Description: Approximate proper motion fitting.
%              This function can not deal with observations crossing the pole. 
% Input  : - Vector of JD per each observation.
%          - Column vector of R.A. [radians] per observation.
%            If Matrix then treat each colum seperatly. 
%          - Column vector of Dec. [radians] per observation.
%            If Matrix then treat each colum seperatly. 
%          - Vector of errors in R.A. [radians] per observation,
%            default is 1.
%            If Matrix then treat each colum seperatly. 
%          - Vector of errors in Dec. [radians] per observation,
%            default is 1.
%            If Matrix then treat each colum seperatly. 
% Output : - Mean epoch for observations.
%          - Best fit parameters:
%            [RA at mean epoch, PM in RA, Dec at mean epoch, PM in Dec]
%            in radians, and radians per days.
%            PM in Declination is corrected for cos(Dec).
%          - Errors in best fit parameters.
%          - \chi^2
%          - Vector of number of degrees of freedom.
%          - Vector of residuals (Observed - Calculated) [radians].
%            The residuals in declination are corrected by cos(Dec).
% Tested : Matlab 7.8
%     By : Eran O. Ofek                   October 2009
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable : 2
%------------------------------------------------------------------------------

if (range(RA)>pi),
   I = find(RA<pi);
   RA(I) = RA(I) + 2.*pi;
end

MeanJD = mean(JD);
N      = length(JD);

Def.ErrRA  = ones(N,1);
Def.ErrDec = ones(N,1);
if (nargin==3),
   ErrRA  = Def.ErrRA;
   ErrDec = Def.ErrDec;
elseif (nargin==4),
   ErrDec = Def.ErrDec;
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end

Nst  = size(RA,2);
Nobs = size(RA,1); 

Hsub = [ones(N,1), JD - MeanJD];
H    = [[Hsub, zeros(N,2)]; [zeros(N,2), Hsub]];

Npar   = 4;
Par    = zeros(Npar,Nst).*NaN;
ParErr = zeros(Npar,Nst).*NaN;
Resid  = zeros(2.*Nobs,Nst).*NaN;
Chi2   = zeros(1,Nst).*NaN;
Dof    = zeros(1,Nst).*NaN;
for Ist=1:1:Nst,
   Y    = [RA(:,Ist);Dec(:,Ist)];
   Ye   = [ErrRA(:,Ist);ErrDec(:,Ist)];

   Ig = find(isnan(Y)==0);  % not NaNs

   if (length(Ig)>3),
      %V      = diag(Ye(Ig).^2);
      %Cov    = inv(H(Ig,:)'*inv(V)*H(Ig,:));
      %Par(:,Ist)    = Cov*H(Ig,:)'*inv(V)*Y(Ig,:);
      %ParErr(:,Ist) = sqrt(diag(Cov));

      [Tmp1,Tmp2] = lscov(H(Ig,:),Y(Ig),1./(Ye(Ig).^2));

      Par(:,Ist) = Tmp1;
      ParErr(:,Ist) = Tmp2;


      Resid(:,Ist) = Y - H*Par(:,Ist);
      Resid(Nobs+1:Nobs.*2,Ist) = Resid(Nobs+1:Nobs.*2,Ist).*cos(Dec(:,Ist));

      Chi2(Ist) = sum((Resid(:,Ist)./Ye).^2);
      Dof(Ist)  = length(Ig) - 4;
   else
      % do nothing
   end
end


% corrected PM for cos(Dec)
Par(4,:) = Par(4,:).*cos(nanmean(Dec));

