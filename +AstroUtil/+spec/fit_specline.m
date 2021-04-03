function [Fit]=fit_specline(Spec,Fun,varargin)
% Fit and measure flux of spectral line
% Package: AstroUtil.spec
% Description: Fit multiple profiles and measure flux of spectral line
% Input  : - Spectrum to fit [Wavelength, Flux].
%            If empty matrix, then will attempt to read the spectrum from
%            the current figure.
%          - Function to fit.
%            If empty matrix than do not fit a function, but calculate
%            line flux and equivalent width.
%          * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'Back'   - Method of background subtractuin:
%                       'none' - use fitting function only.
%                       'man'  - Interactive subtraction of background
%                                (default).
%                       Alternatively, this can be a vector that specify
%                       from where to obtain the background:
%                       [X1, X2, X3, X4].
%                       Where X1 and X2 mark the position of the left
%                       background region and X3 and X4 is for the right
%                       background region.
%            'Poly'   - Polynomial order of manual background subtraction.
%                       Default is 1.
%            'Plot'   - Plot fitted lines {'y'|'n'}. Default is 'y'.
%            'PlotPar'- Cell array of additional parameters to pass to the
%                       fitted lines ploting function.
%                       Default is {'--','Color',[0.8 0.8 0.8]}.
%            'Par0'   - Vector of guess parameters.
%                       If equal 'guess' then will attempt to call the
%                       fitted function Fun('guess',Spec) which suppose
%                       to return a guess.
%            'Nsim'   - Number of simulations within the specified
%                       background ranges. Default is 1.
%            'MinBkg' - Minimum bkg range in terms of fraction of the defined
%                       range. Default is 0.3 (30%).
%
% Output :   - Vector of Fit structure.
%              If Nsim=1 then the first Fit values (Flux,EW) are according to the original
%              background ranges. The rest are according to the different
%              randomly-selected bkg ranges, but always the range from
%              which the flux is calculated (Iall) is of the inner edges of
%              the original bkg (supplied or by interactive selection).
% Tested : Matlab 7.13
%     By : Eran Ofek / Ofer Yaron          Jul 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Fit=fit_specline(Spec,@fun_gauss);
%          Fun = @(Par,X) fun_gauss(Par(1:3),X)+fun_gauss(Par(4:6),X)
%          Fit=fit_specline(Spec,Fun,'Par0',[1 6630 10 1 6630 100]);
%          Fun = @(Par,X) fun_gauss(Par(1:3),X)+fun_lorentzian(Par(4:6),X)
%          Fit=AstroUtil.spec.fit_specline(Spec,Fun,'Par0',[1 6630 10 6630 100 10]);
%          Fit=fit_specline(Spec,[],'Back',bkg,'Nsim',100,'MinBkg',0.5);
% Reliable: 2
%------------------------------------------------------------------------------

import Util.fit.*

DefV.Back    = 'man';
DefV.Poly    = 1;
DefV.Plot    = 'y';
DefV.PlotPar = {'--','Color',[0.8 0.8 0.8]};
DefV.Par0    = 'guess';
DefV.Nsim    = 1;
DefV.MinBkg  = 0.3;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if nargin<2
    Fun = [];
end

OpenNewPlot = true;
if (isempty(Spec))
   Spec = [get(get(gca,'Children'),'XData'), get(get(gca,'Children'),'YData')];
   OpenNewPlot = false;
end

% DBG:
%clrmat=jet;
%clridx = 1:floor(length(clrmat)/InPar.Nsim):length(clrmat);

for j=1:InPar.Nsim

    if (~ischar(InPar.Back))
       % assume InPar.Back is a 4 elements vector
       if j==1  % First bkg ranges are always the default/supplied ones
          XL = InPar.Back(1:2);
          XR = InPar.Back(3:4);
          XLorig=XL;
          XRorig=XR;
       else
           [XL XR]=set_random_bkg(XLorig, XRorig, InPar.MinBkg);
       end
       IL  = find(Spec(:,1)>min(XL) & Spec(:,1)<max(XL));
       IR  = find(Spec(:,1)>min(XR) & Spec(:,1)<max(XR));
       %Iall= find(Spec(:,1)>(min( min(XL),min(XR) )) & ...
       %           Spec(:,1)<(max( max(XL),max(XR) )) );
       Iall= find(Spec(:,1)>(min( max(XLorig),max(XRorig) )) & ...
                  Spec(:,1)<(max( min(XLorig),min(XRorig) )) );
       %The above is similar to the line below just more error proof - left vs right
       %Iall= find(Spec(:,1)>=max(XLorig) & Spec(:,1)<=min(XRorig));       
       Ib  = [IL;IR];
       Res = fitgenpoly(Spec(Ib,1),Spec(Ib,2),1,InPar.Poly);   
       Back = [Spec(Iall,1),polyval(Res.Par,Spec(Iall,1))];

    else
       switch lower(InPar.Back)
        case 'man'
            if j==1
               % interactive background subtraction
               if (OpenNewPlot)
                  stairs(Spec(:,1),Spec(:,2),'k-');
               end

               ReFit = true;
               Hpb   = [];
               while (ReFit)
                  fprintf('Zoom on\n');
                  zoom on;
                  input('press any key to continue','s');
                  zoom off;
                  fprintf('Select background region left of line to fit\n');
                  fprintf('Two mouse left click to select background region\n');
                  [XL,YL] = ginput(2);
                  fprintf('Select background region right of line to fit\n');
                  fprintf('Two mouse left click to select background region\n');
                  [XR,YR] = ginput(2);
                  
                  XLorig=XL;
                  XRorig=XR;

                  % calc background
                  IL  = find(Spec(:,1)>min(XL) & Spec(:,1)<max(XL));
                  IR  = find(Spec(:,1)>min(XR) & Spec(:,1)<max(XR));
                  %Iall= find(Spec(:,1)>(min( min(XL),min(XR) )) & ...
                  %           Spec(:,1)<(max( max(XL),max(XR) )) );
                  Iall= find(Spec(:,1)>(min( max(XLorig),max(XRorig) )) & ...
                             Spec(:,1)<(max( min(XLorig),min(XRorig) )) );
                  %The above is similar to the line below just more error proof - left vs right
                  %Iall= find(Spec(:,1)>=max(XLorig) & Spec(:,1)<=min(XRorig));       
                  Ib  = [IL;IR];
                  Res = fitgenpoly(Spec(Ib,1),Spec(Ib,2),1,InPar.Poly);
                  if (~isempty(Hpb))
                     delete(Hpb);
                  end
                  Back = [Spec(Iall,1),polyval(Res.Par,Spec(Iall,1))];
                  switch lower(InPar.Plot)
                   case 'y'
                      hold on;
                      Hpb = plot(Back(:,1),Back(:,2),InPar.PlotPar{:});
                   otherwise
                      % do nothing
                  end

                  R = input('Do you want to re-select/fit background (Y/[N]) : ','s');
                  switch lower(R)
                   case 'y'
                        ReFit = true;
                   otherwise
                        ReFit = false;
                  end
               end
            else
               [XL XR]=set_random_bkg(XLorig, XRorig, InPar.MinBkg);

               IL  = find(Spec(:,1)>min(XL) & Spec(:,1)<max(XL));
               IR  = find(Spec(:,1)>min(XR) & Spec(:,1)<max(XR));
               %Iall= find(Spec(:,1)>(min( min(XL),min(XR) )) & ...
               %           Spec(:,1)<(max( max(XL),max(XR) )) );
               Iall= find(Spec(:,1)>(min( max(XLorig),max(XRorig) )) & ...
                          Spec(:,1)<(max( min(XLorig),min(XRorig) )) );
               %The above is similar to the line below just more error proof - left vs right
               %Iall= find(Spec(:,1)>=max(XLorig) & Spec(:,1)<=min(XRorig));       
               Ib  = [IL;IR];
               Res = fitgenpoly(Spec(Ib,1),Spec(Ib,2),1,InPar.Poly);
               Back = [Spec(Iall,1),polyval(Res.Par,Spec(Iall,1))];
            end
        case 'none'
           % do nothing
           Back = [Spec(:,1), zeros(size(Spec(:,1)))];
           Iall = (1:1:size(Spec,1)).';
        otherwise
           error('Unknown Back option');
       end
    end

    % DBG:
    %clr=clrmat(clridx(j),:);
    %plot(Back(:,1),Back(:,2),'color',clr);    

    if (~isempty(Fun))
       % fit
       if (~ischar(InPar.Par0))
          Par0 = InPar.Par0;
       else
          try
             Par0 = feval(Fun,'guess',Spec(Iall,:));
          catch
             error('Guess parameters (Par0) is required');
          end
       end

       [Fit.Par,Fit.Resid,Fit.J,Fot.Cov,Fit.Err] = nlinfit_my(Spec(Iall,1),Spec(Iall,2)-Back(:,2),Fun,Par0);
       ModelY = feval(Fun,Fit.Par,Spec(Iall,1));
       switch lower(InPar.Plot)
        case 'y'
           hold on;
           plot(Spec(Iall,1),ModelY+Back(:,2),InPar.PlotPar{:});
       otherwise
           % do nothing
       end
       % calculate lines parameters
       Fit(j).ModelTotalFlux = trapz(Spec(Iall,1),ModelY);
       Fit(j).ModelTotalEW   = trapz(Spec(Iall,1),1 - ModelY./Back(:,2));
       Fit(j).ObsTotalFlux   = trapz(Spec(Iall,1),Spec(Iall,2)-Back(:,2));
       Fit(j).ObsTotalEW     = trapz(Spec(Iall,1),1 - Spec(Iall,2)./Back(:,2));
       Fit(j).Back           = Back;
       Fit(j).ModelY         = ModelY;
    else
       % do not fit  calcl Flux and EW
       Fit(j).ObsTotalFlux   = trapz(Spec(Iall,1),Spec(Iall,2)-Back(:,2));
       Fit(j).ObsTotalEW     = trapz(Spec(Iall,1),1 - Spec(Iall,2)./Back(:,2));
       Fit(j).Back           = Back;

    end



    if (1==0)

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

        if (length(Par0)>4)
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
    if (length(Par0)>3)
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
    
end

end

%------------------------------------------------------
function [XL XR]=set_random_bkg(XLorig, XRorig, minbkg)
% choose random background

   LRange=XLorig(2)-XLorig(1); % Right Range
   RRange=XRorig(2)-XRorig(1); % Left Range

   ReFit = true;
   % Looping till the new (randomly obtained) ranges 
   % match the required MinBkg criterion.
   while (ReFit)
      frac=[rand rand];
      frac=sort(frac);
      XL=[XLorig(1)+frac(1)*LRange XLorig(1)+frac(2)*LRange];
      XR=[XRorig(1)+frac(1)*RRange XRorig(1)+frac(2)*RRange];
      XL=round(XL);
      XR=round(XR);
      LRange_new=XL(2)-XL(1);
      RRange_new=XR(2)-XR(1);
      if (LRange_new>=minbkg*LRange && RRange_new>=minbkg*RRange)
          ReFit = false;
      end
   end
end
%------------------------------------------------------

end
