function [RMS,Dist,ImPos,ConvInfo]=find_images_regions(BetaPos,ModelPars,ModelType,StartImPos,Dls_Ds,MaxDist,SearchMethod,ConvThreshold);
%-----------------------------------------------------------------------------
% find_images_regions function                                          glens
% Description:  Given a mass model, its parameters and a
%                               source position, search for images of
%                               the source only in a predefind given regions
%                               or starting points.
% Input  : - A single source positions [X, Y].
%          - Model parameters (see: calc_alpha.m).
%          - Model type (see: calc_alpha.m).
%          - A single source positions [X, Y].
%          - "Guess" images position to start the search in, [X, Y]
%            or [X, Y, ErrX, ErrY] (in that case return chi2 instead of RMS).
%          - Dls/Ds
%          - Maximum distance to search within -
%            default is NaN.
%          - SearchMethod:
%            'Jacobian'   - Use the Jacobian to converge, default.
%            'Direct'     - Use 'fminsearch' (local minimum)
%            'Scan'       - Scan a region with a resolution given by Resolution.
%          - Convergence threshold for solution in actual units,
%            For 'Jacobian' default is 0.01.
%            For 'Scan'     default is 1.
% Output : - RMS or Chi2.
%            Return NaN if solution not found within MaxNoIter iterations.
%          - Vector of distances between the observed images position
%            and the calculated images position
%            (x components followed by y components).
%          - Images position and Jacobian, corresponding to the source
%            [X, Y, A_11, A_22, A_12].
%            If an images wasn't found within the search distance
%            (or not converged) set the image position to NaN.
%          - Vector of convergence: 0 - if not converged; 1 - if converged
% Tested : Matlab 7.0
%     By : Eran O. Ofek        April 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------------------
MaxNoIter = 1;
ColNorm   = 7;

if (nargin==5),
   MaxDist       = NaN;
   SearchMethod  = 'Jacobian';
   ConvThreshold = [];
elseif (nargin==6),
   SearchMethod  = 'Jacobian';
   ConvThreshold = [];
elseif (nargin==7),
   ConvThreshold = [];
elseif (nargin==8),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(ConvThreshold)==1),
   switch SearchMethod
    case 'Jacobian'
       ConvThreshold = 0.01;
    case 'Scan'
       ConvThreshold = 1.0;
    case 'Direct'
       ConvThreshold = [];
    otherwise
       error('Unknown SearchMethod Option');
   end
end

Nim = size(StartImPos,1);

switch SearchMethod
 case 'Jacobian'
    %---
    % Use the Jacobian terms A_11, A_22, A_12 to
    % solve linearly for Theta(Beta)
    %---

    ImPos        = zeros(Nim,5);
    ImPos(:,1:2) = StartImPos(:,1:2);
    Theta_Beta   = zeros(Nim,2);
    Diff         = zeros(Nim,2);
    %ConvInfo     = zeros(Nim,1);
    ConvInfo     = 1;

    Converge = 0;
    Iconv    = 0;
    while (Converge==0),
       Iconv = Iconv + 1;
       MP_DlsDs = ModelPars;
       MP_DlsDs(:,ColNorm) = MP_DlsDs(:,ColNorm).*Dls_Ds;
       [AlphaX,AlphaY,A_11,A_22,A_12] = calc_alpha(ImPos(:,1),ImPos(:,2),MP_DlsDs,ModelType);
       %AlphaX = AlphaX.*Dls_Ds;
       %AlphaY = AlphaY.*Dls_Ds;

       %----------------------
       %--- For each image ---
       %----------------------
       %--- X solution ---
       % bx0 -> BetaGridMatX
       % by0 -> BetaGridMatY
       % ax0 -> AlphaX(Iim)
       % ay0 -> AlphaY(Iim)
       % x0  -> ThetaX(Iim)
       % y0  -> ThetaY(Iim)
       % -(-by0*Dxy-bx0-ax0+Dxx*x0+Dxy*y0-ay0*Dxy+Dyy*bx0+Dyy*ax0-Dyy*Dxx*x0+Dyx*x0*Dxy)/
       %  (1-Dxx-Dyy+Dyy*Dxx-Dyx*Dxy)
       %--- Y solution ---
       %  ( bx0*Dyx-Dxy*y0*Dyx-Dxx*by0-Dyx*x0+by0+ay0+ax0*Dyx-Dyy*y0-Dxx*ay0+y0*Dyy*Dxx)/
       %  (1-Dxx-Dyy+Dyy*Dxx-Dyx*Dxy)
   
       Det = 1 - A_11 - A_22 + A_22.*A_11 - A_12.^2;
   
       %  ThetaX(Beta)
       Theta_Beta(:,1) = -(-BetaPos(2) .* A_12 - BetaPos(1) - AlphaX + A_11 .* ImPos(:,1) + A_12 .* ImPos(:,2) -  AlphaY .* A_12 + A_22 .* BetaPos(1) + A_22 .* AlphaX - A_22.*A_11 .* ImPos(:,1) + A_12 .* ImPos(:,1) .* A_12)./Det;
   
       %  ThetaY(Beta)
       Theta_Beta(:,2) = (BetaPos(1) .* A_12 - A_12 .* ImPos(:,2) .* A_12 - A_11 .* BetaPos(2) - A_12 .* ImPos(:,1) + BetaPos(2) + AlphaY + AlphaX .* A_12 - A_22 .* ImPos(:,2) - A_11 .* AlphaY + ImPos(:,2) .* A_22 .* A_11)./Det;
   
       Diff         = Theta_Beta - ImPos(:,1:2);
       ImPos(:,1:2) = Theta_Beta;



       if (max(sqrt(sum(Diff.^2,2)))<ConvThreshold | Iconv>=MaxNoIter),
          Converge = 1;
       end

    end
    if (Iconv>MaxNoIter),
       %--- Not converged ---
       %ImPos  = NaN.*ImPos;
       ConvInfo = 0;
    else
       ImPos(:,3:5) = [A_11, A_22, A_12];
    end

 case 'Direct'
    error('Direct method is not available');
 case 'Scan'
    %---
    % Search the images position by calculating the deflection for a region
    % around each guess image position...
    %---

    ImPos = zeros(Nim,5);
    [GridDiffX, GridDiffY] = meshgrid([-MaxDist:ConvThreshold:MaxDist]);
    Size = size(GridDiffX);

    for Iim=1:1:Nim,

       MP_DlsDs = ModelPars;
       MP_DlsDs(:,ColNorm) = MP_DlsDs(:,ColNorm).*Dls_Ds;
       [AlphaX,AlphaY,A_11,A_22,A_12] = calc_alpha(StartImPos(Iim,1)+GridDiffX,StartImPos(Iim,2)+GridDiffY,MP_DlsDs,ModelType);
       %AlphaX = AlphaX.*Dls_Ds;
       %AlphaY = AlphaY.*Dls_Ds;
       %--- calculate Beta = Theta - Alpha(Theta) ---
       BetaX      = StartImPos(Iim,1) + GridDiffX - AlphaX;
       BetaY      = StartImPos(Iim,2) + GridDiffY - AlphaY;
       %--- Compare Calculated Beta with given Beta (source position) ---
       %surface([-MaxDist:ConvThreshold:MaxDist],[-MaxDist:ConvThreshold:MaxDist],sqrt((BetaPos(1)-BetaX).^2 + (BetaPos(2)-BetaY).^2));shading interp; colorbar

       [Min,MinI] = min2d(sqrt((BetaPos(1)-BetaX).^2 + (BetaPos(2)-BetaY).^2));
       ImPos(Iim,1) = StartImPos(Iim,1) + GridDiffX(MinI(1),MinI(2));
       ImPos(Iim,2) = StartImPos(Iim,2) + GridDiffY(MinI(1),MinI(2));
       if (MinI(1)==1 | MinI(2)==1 | MinI(1)==Size(1) | MinI(2)==Size(1)),
          %--- Solution isoutside search region ---
          ImPos = ImPos.*NaN;
       else
          ImPos(Iim,3) = A_11(MinI(1),MinI(2));
          ImPos(Iim,4) = A_22(MinI(1),MinI(2));
          ImPos(Iim,5) = A_12(MinI(1),MinI(2));
       end

    end
 otherwise
    error('Unknown SearchMethod Option');
end


Dist = [StartImPos(:,1) - ImPos(:,1); StartImPos(:,2) - ImPos(:,2)];
if (size(StartImPos,2)>2),
   %--- X/Y errors given ---
   %--- claculate chi2   ---
   RMS = Dist./[StartImPos(:,3); StartImPos(:,4)];
else
   %--- No X/Y errors ---
   %--- calculate RMS ---
   RMS  = sqrt(mean(Dist.^2));
end

