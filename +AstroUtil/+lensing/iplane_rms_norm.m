function [RMS_I,AllDist_I,Dof,AllBestBeta]=iplane_rms_norm(Norm,ModelPars,ModelType,ImagesCell,Dls_Ds);
%--------------------------------------------------------------------------
% iplane_rms_norm function                                           glens
% Description:   Given a model parameters of a lens, and the
%                          images position, calculate the best source
%                          position that minimize the residuals in
%                          image plane (see also: iplane_rms_smart.m).
%                          ================================================
%                          This program is a duplicate of iplane_rms.m
%                          The only difference is that this program get
%                          the lens model normalization (for the first line
%                          in the model - the rest of the lines are
%                          normalized proportianly) as the first
%                          input argument.
%                          ================================================
% Input  : - Model Normalization (override the 7th column in ModelPars). 
%          - Model parameters (see calc_alpha.m).
%          - Vector of Model Type (see calc_alpha.m).
%          - Images position cell array:
%            Each cell contains [ThetaX, ThetaY, ErrorX, ErrorY]
%            of the images of a single source [pixels units!].
%            It is recomended that the first image in each list
%            will be the faintest image corresponds to each source.
%          - Vector of Dls/Ds for each source.
% Output : - RMS or chi2 for best image plane solution.
%            In case image plane solution not attempted set to NaN.
%          - Residuals in the image plane (x components for all images
%            of a source followed by the y component, followed
%            by next source...)
%          - Number of constraints.
%          - Best source position [X, Y].
%            In case the program didn't solve for source position
%            return [NaN, NaN].
% Tested : Matlab 7.0
%     By : Eran Ofek              April 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See Also: iplane_rms.m, iplane_rms_smart.m
%--------------------------------------------------------------------------

GuessBetaMethod        = 'MagIsoWeighted';  % {'FirstImage' | 'Mean' | 'Median' | 'MagIsoWeighted'}

ColNorm     = 7;

%--- noramlized Norm according to first line ---
NormalizedNorm        = ModelPars(:,ColNorm)./ModelPars(1,ColNorm);
%--- Set the Normalizations for all lines ---
ModelPars(:,ColNorm)  = Norm .* NormalizedNorm;



%------------------------------
%--- For each set of images ---
%------------------------------
AllDist_I  = zeros(0,0);
AllErr_I   = zeros(0,0);
AllBestBeta= zeros(0,2);
Nsource    = length(ImagesCell);
for Is=1:1:Nsource,
   ObsImPos = ImagesCell{Is};
   Nim      = size(ObsImPos,1);
   %------------------------------------------------------------------
   %--- Calculate deflection and magnification for images position ---
   %------------------------------------------------------------------
   [AlphaX,AlphaY,A_11,A_22,A_12] = calc_alpha(ObsImPos(:,1),ObsImPos(:,2),ModelPars,ModelType);
   %--- Calculate source position corresponding to each image position ---
   BetaAllIm = [ObsImPos(:,1) - Dls_Ds(Is).*AlphaX, ObsImPos(:,2) - Dls_Ds(Is).*AlphaY];
      

   %----------------------------------------------------
   %--- Calculate approximate source position (Beta) ---
   %----------------------------------------------------
   %--- Method to use in calculating approximate source position ---
   switch GuessBetaMethod
    case 'FirstImage'
       %--- Use first image to estimat source position (Beta) ---
       % recomanded to use the smallest magnification image
       [AlphaX,AlphaY,A_11,A_22,A_12] = calc_alpha(ObsImPos(1,1),ObsImPos(1,2),ModelPars,ModelType);
       BetaStart = BetaAllIm(1,:);
      
    case 'Mean'
       %--- Use all images to estimate source positions and use the mean (Beta) ---
       BetaStart = mean(BetaAllIm);
      
    case 'Median'
       %--- Use all images to estimate source positions and use the mean (Beta) ---
       BetaStart = median(BetaAllIm);
      
    case 'MagIsoWeighted'
       %--- Weight by magnitude (isotropic) ---
       Mag = 1./((1 - A_11).*(1 - A_22) - A_12.^2);
       BetaStart = [wmean([BetaAllIm(:,1), 1./Mag]), wmean([BetaAllIm(:,2), 1./Mag])];
      
    otherwise
       error('Unknown GuessBetaMethod');
   end
      
      
   %---------------------------------------------------------------------------------------
   %--- Solve for sorce position (Beta) that is minimize the scatter in the image plane ---
   %---------------------------------------------------------------------------------------
   MinimizeMethod = 'Direct';
   [BestBeta,RMS_I,Dist_I,ImPos]=beta_minimize(BetaStart,ModelPars,ModelType,ObsImPos,Dls_Ds(Is),MinimizeMethod);   
   AllBestBeta= [AllBestBeta; BestBeta];
   AllDist_I  = [AllDist_I;  Dist_I];
   % all corresponding image position errors
   if (size(ObsImPos,2)>2),
      AllErr_I   = [AllErr_I;   ObsImPos(:,3); ObsImPos(:,4)];
   else
      AllErr_I   = NaN;
   end
end


%--------------------------
%--- Calculate RMS/Chi2 ---
%--------------------------
if (size(ObsImPos,2)>2),
   %--- Calculate chi2 ---
   RMS_I = sum( (AllDist_I./AllErr_I).^2 );
   Dof   = length(AllDist_I);
else
   %--- Calculate RMS ---
   RMS_I = std(AllDist_I);
   Dof   = length(AllDist_I);
end


