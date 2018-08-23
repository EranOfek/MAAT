function [Fit,FittedSpec]=fit_specline(Spec,Fun);
%------------------------------------------------------------------------------
% fit_specline function                                              AstroSpec
% Description:
% Input  : -
% Output : -
% Tested : Matlab 7.13
%     By : Eran Ofek                         July 2012
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%------------------------------------------------------------------------------

Def.Interactive = 'y';
%if (nargin==1),
   Interactive = Def.Interactive;
%end

switch lower(Interactive)
 case 'y'
    stairs(Spec(:,1),Spec(:,2),'k','LineWidth',1)

    Norm = nanmedian(Spec(:,2));
    Spec(:,2) = Spec(:,2)./Norm;
    Fit.Norm = Norm;

    % take off NaNs
    I = find(~isnan(Spec(:,2)));
    Spec = Spec(I,:);


%    stairs(Spec(:,1),Spec(:,2),'k','Color',[0.8 0.8 0.8],'LineWidth',2)
    hold on;

    zoom on;
    R = input('Zoom in on region to fit - press return to continue','s');

    fprintf('Mark continuum on left side of line (using two mouse clicks)\n');
    [ContLX,~] = ginput(2);
    Icl = find(Spec(:,1)>=min(ContLX) & Spec(:,1)<=max(ContLX));
    fprintf('Mark continuum on right side of line (using two mouse clicks)\n');
    [ContRX,~] = ginput(2);
    Icr = find(Spec(:,1)>=min(ContRX) & Spec(:,1)<=max(ContRX));
    Ic  = [Icl;Icr];
    Irange = [min(Ic):1:max(Ic)];
   
Deg=1
    Par = fitpoly(Spec(Ic,1),Spec(Ic,2),ones(size(Ic)),Deg,3);
    ContY = polyval(flipud(Par),Spec(Irange,1));
    plot(Spec(Irange,1),ContY,'b--');

    fprintf('Best fit contimuum is shown in blue line\n');


%fit_gaussians
%    Fun = @(Par,WL) Par(1).*exp(-(WL-Par(2)).^2./(2.*Par(3).^2))
%    Fun = @(Par,WL) Par(1).*exp(-(WL-Par(2)).^2./(2.*Par(3).^2)) + Par(4).*exp(-(WL-Par(5)).^2./(2.*Par(6).^2))
%     Fun = @(Par,WL) Par(1).*exp(-(WL-Par(2)).^2./(2.*Par(3).^2)) + Par(4).*exp(-(WL-Par(5)).^2./(2.*Par(6).^2)) + Par(7).*exp(-(WL-Par(8)).^2./(2.*Par(9).^2))

    % get best guess
    Par0    = feval(Fun,'get_guess',Spec(Irange,:));
    
%    Par0(1) = mean(Spec(Irange,2));  % amplitude
%    Par0(2) = mean(Spec(Irange,1));  % central wavelength
%    Par0(3) = (max(Spec(Irange,1)) - min(Spec(Irange,1)))./10; % std
%    Par0(4) = Par0(1).*2;
%    Par0(5) = Par0(2);
%    Par0(6) = Par0(3)./4;
%    Par0(7) = Par0(1).*3;
%    Par0(8) = Par0(2);
%    Par0(9) = Par0(3)./8;


    WL      = Spec(Irange,1);
    %BestPar = nlinfit(WL, Spec(Irange,2)-ContY,Fun,Par0);
    %BestPar = nlinfit(WL, Spec(Irange,2)-ContY,@fun_2gauss_conv,Par0);
    [BestPar,Fit.Resid,~,Fit.Cov] = nlinfit(WL, Spec(Irange,2)-ContY,Fun,Par0);
    Fit.TotFlux  = trapz(WL,feval(Fun,BestPar,WL) + ContY);
    Fit.ContFlux = trapz(WL,ContY);
    
    if (length(Par0)>4),
       Fit.Flux1    = trapz(WL, feval(Fun,[BestPar(1:3), zeros(1,3)],WL));
       Fit.Flux2    = trapz(WL, feval(Fun,[zeros(1,3), BestPar(4:6)],WL));
    else
       Fit.Flux1    = NaN;
       Fit.Flux2    = NaN;
    end
    Fit.EqW      = trapz(WL, 1 - (feval(Fun,BestPar,WL) + ContY)./ContY);

    Fit.BestPar = BestPar;   


% renormalize
Fit.BestPar(1) = Fit.BestPar(1).*Fit.Norm;
if (length(Par0)>3),
   Fit.BestPar(4) = Fit.BestPar(4).*Fit.Norm;
end
Fit.TotFlux    = Fit.TotFlux.*Fit.Norm;
Fit.ContFlux   = Fit.ContFlux.*Fit.Norm;
Fit.Flux1      = Fit.Flux1.*Fit.Norm;
Fit.Flux2      = Fit.Flux2.*Fit.Norm;
Fit.StdResid   = std(Fit.Resid).*Fit.Norm;

FittedSpec = [WL,Norm.*(feval(Fun,BestPar,WL) + ContY)];
plot(WL,Norm.*(feval(Fun,BestPar,WL) + ContY),'g--','Color',[0.8 0.8 0.8],'LineWidth',2);




%SubContSpec = [Spec(:,1), Spec(:,2)-ContY




 otherwise


end
