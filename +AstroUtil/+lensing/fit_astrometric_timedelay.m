function [Res]=fit_astrometric_timedelay(T,F,X,X1,X2,varargin)
% SHORT DESCRIPTION HERE
% Package: AstroUtil.lensing
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res]=AstroUtil.lensing.fit_astrometric_timedelay
% Reliable: 
%--------------------------------------------------------------------------

Def.X1 = 0;
Def.X2 = 1;
SimMode = false;
if (nargin==0)
    SimMode = true;
    X1 = Def.X1;
    X2 = Def.X2;
elseif (nargin==3)
    X1 = Def.X1;
    X2 = Def.X2;
else
    % do nothing
end

if (~SimMode)
    if (isempty(F))
        SimMode = true;
    end
end


DefV.TauVec               = (0:1:200)'; %logspace(0,log10(200),200).'; %(0:1:200)';     % trial Tau (time delay)
DefV.AlphaVec             = (0.04:0.02:1)';   % trial Alpha =f2(t)/f1(t+Tau)
DefV.InterpMethod         = 'linear';
DefV.D_X1                 = 0.001;
DefV.D_X2                 = 0.001;
DefV.D_F                  = 0.03;
DefV.D_X                  = 0.03;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (SimMode)
    % simulation mode
    Alpha0      = 0.2;   % the real Alpha
    Tau0        = 50;    % the real Tau
    MaxTime     = 500;
    Ta          = (1:1:MaxTime+Tau0)';
    %X1          = 0;
    %X2          = 1;
    PowerLawInd = 2;
    Norm        = -0.3;   % std of time series
    Err         = 0.03;
    % power-law time series (in magnitude)
    f  = Util.stat.rand_ps(Ta,[PowerLawInd, Norm]);
    f  = f(:,2);   % value
    MeanMag = 18;
    Mag  = f + MeanMag + randn(size(f)).*Err;
    % convert to flux
    ZP = 22;
    f  = 10.^(-0.4.*(Mag-ZP));
    f1 = f(1:MaxTime);
    f2 = Alpha0.*f(1+Tau0:end);  % = f1(t+Tau)
    T  = Ta(1:MaxTime);
    F  = f1+f2;
    ErrPos = 0.03;
    X  = (X1.*f1 + X2.*f2)./F + randn(size(f1)).*ErrPos;
    plot(T,f1);
    hold on;
    plot(T,f2);
    plot(T,X);
    
end




D_X1 = InPar.D_X1;
D_X2 = InPar.D_X2;
D_F  = InPar.D_F;
D_X  = InPar.D_X;

Ntau   = numel(InPar.TauVec);
Nalpha = numel(InPar.AlphaVec);

Res.RMS      = zeros(Nalpha,Ntau);
Res.TauVec   = InPar.TauVec;
Res.AlphaVec = InPar.AlphaVec;
Res.Chi2     = zeros(Nalpha,Ntau);
Res.Dof      = zeros(Nalpha,Ntau);
Res.FluxRMS  = nanstd(F);
Res.AstrometricRMS = nanstd(X);
Res.Npar     = 2;
for Ialpha=1:1:Nalpha
    Alpha = InPar.AlphaVec(Ialpha);
    % estimate f1(t+Tau)
    %f = (X.*F - X1.*F)./(Alpha.*(X2-X1));
    f = F.*(X - X2)./(X1 - X2);
    %[ErrExp,ErrVar]=Util.symbolic.symerror((X*F - X2*F)/(Alpha*(X1-X2)),X,F,X1,X2)
    %D_f = (D_X1.^2.*(F./(Alpha.*(X1 - X2)) + (F.*X - F.*X1)./(Alpha.*(X1 - X2).^2)).^2 + (D_X.^2.*F.^2)./(Alpha.^2.*(X1 - X2).^2) + (D_X2.^2.*(F.*X - F.*X1).^2)./(Alpha.^2.*(X1 - X2).^4) + (D_F.^2.*(X - X1).^2)./(Alpha.^2.*(X1 - X2).^2)).^(1./2);
    
    f = F.*(X - X2)./(X1 - X2);
    %[ErrExp,ErrVar]=Util.symbolic.symerror((X*F - X2*F)/(X1-X2),X,F,X1,X2)
    D_f = (D_X2.^2.*(F./(X1 - X2) - (F.*X - F.*X2)./(X1 - X2).^2).^2 + (D_X.^2.*F.^2)./(X1 - X2).^2 + (D_X1.^2.*(F.*X - F.*X2).^2)./(X1 - X2).^4 + (D_F.^2.*(X - X2).^2)./(X1 - X2).^2).^(1./2);
    
    for Itau=1:1:Ntau
        % for each trial time delay
        Tau      = InPar.TauVec(Itau);
        f1       = interp1(T-Tau,f,T,InPar.InterpMethod);
        if (numel(D_f)==1)
            D_f1 = D_f;
        else
            D_f1     = interp1(T+Tau,D_f,T,'nearest'); 
        end
        %Fmodel   = f1 + Alpha.*f;
        Fmodel   = f + Alpha.*f1;
        
        D_Fmodel = (D_f1.^2 + Alpha.^2.*D_f.^2).^(1./2);
        
        DeltaFF               = F - Fmodel;
        Res.RMS(Ialpha,Itau)  = nanstd(DeltaFF);
        Res.Chi2(Ialpha,Itau) = nansum((DeltaFF).^2./(D_F.^2 + D_Fmodel.^2));
        Res.Dof(Ialpha,Itau)  = sum(~isnan(DeltaFF)) - Res.Npar;
    end
end

[MinVal,MinInd] = Util.stat.minnd(Res.Chi2);
Res.BestAlpha   = Res.AlphaVec(MinInd(1));
Res.BestTau     = Res.TauVec(MinInd(2));
Res.BestChi2    = Res.Chi2(MinInd(1),MinInd(2));
Res.BestDof     = Res.Dof(MinInd(1),MinInd(2));

% surface(Res.TauVec,Res.AlphaVec,Res.Chi2); colorbar
%contour(Res.TauVec,Res.AlphaVec,Res.Chi2,min(Res.Chi2(:))+chi2inv([0.68 0.95 0.9973],2)); colorbar
