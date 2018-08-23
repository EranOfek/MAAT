function [RMS_I,RMS_PI,RMS_S,AllDist_I,AllDist_PI,AllDist_S,Dof,AllBestBeta]=iplane_rms_smart(ModelPars,ModelType,ImagesCell,Dls_Ds);
%--------------------------------------------------------------------------
% iplane_rms_smart function                                          glens
% Description: Given a model parameters of a lens, and the
%                          images position, calculate the best source
%                          position that minimize the residuals in
%                          image plane.
%                          The program starts with source plane fitting
%                          and pseudo source plane fitting - if
%                          successful proceed with image plane fitting.
% Input  : - Model parameters (see calc_alpha.m).
%          - Vector of Model Type (see calc_alpha.m).
%          - Images position cell array:
%            Each cell contains [ThetaX, ThetaY, ErrorX, ErrorY]
%            of the images of a single source [pixels units!].
%            It is recomended that the first image in each list
%            will be the faintest image corresponds to each source.
%          - Vector of Dls/Ds for each source.
% Output : - RMS or chi2 for best image plane solution.
%            In case image plane solution not attempted set to NaN.
%          - RMS for pseudo image plane solution.
%            In case image plane solution not attempted set to NaN.
%          - RMS for source plane solution.
%          - Residuals in the image plane (x components for all images
%            of a source followed by the y component, followed
%            by next source...).
%          - Residuals in the pseudo image plane (x components for all images
%            of a source followed by the y component, followed
%            by next source...).
%            In case image plane solution was achived return NaN.
%          - Residuals in the source plane (x components for all images
%            of a source followed by the y component, followed
%            by next source...).
%            In case image plane solution was achived return NaN.
%          - Number of constraints.
%          - Best source position [X, Y].
%            In case the program didn't solve for source position
%            return [NaN, NaN].
% Tested : Matlab 7.0
%     By : Eran Ofek              April 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See Also: iplane_rms.m, iplane_rms_norm.m
%--------------------------------------------------------------------------

CheckSourcePosScat     = 'StD';             % {'StD' | 'MaxMean' | 'MaxMedian' | 'no'}
CheckSingleBeta        = 'yes';             % {'yes' | 'no'}
SourcePosScatCrit      = 40;                % Source Plane Scatter <--- asuming pixel units
PseImagePosScatCrit    = 40;                % Pseudo Image Plane Scatter <--- asuming pixel units
GuessBetaMethod        = 'MagIsoWeighted';  % {'FirstImage' | 'Mean' | 'Median' | 'MagIsoWeighted'}
StopIf_S_catterLarge   = 0;                 % {0 | 1}
StopIf_PI_ScatterLarge = 0;                 % {0 | 1}
StopIf_I_ScatterLarge  = 0;                 % {0 | 1}

%------------------------------
%--- For each set of images ---
%------------------------------
AllDist_S  = zeros(0,0);
AllDist_PI = zeros(0,0);
AllDist_I  = zeros(0,0);
AllErr_I   = zeros(0,0);
AllBestBeta= zeros(0,2);
Nsource    = length(ImagesCell);
Stop       = 0;
for Is=1:1:Nsource,
   if (Stop==0),
      ObsImPos = ImagesCell{Is};
      Nim      = size(ObsImPos,1);
      %------------------------------------------------------------------
      %--- Calculate deflection and magnification for images position ---
      %------------------------------------------------------------------
      [AlphaX,AlphaY,A_11,A_22,A_12] = calc_alpha(ObsImPos(:,1),ObsImPos(:,2),ModelPars,ModelType);
      %--- Calculate source position corresponding to each image position ---
      BetaAllIm = [ObsImPos(:,1) - Dls_Ds(Is).*AlphaX, ObsImPos(:,2) - Dls_Ds(Is).*AlphaY];
      
      %----------------------------------------
      %--- test for scatter in source plane ---
      %----------------------------------------
      % if scatter in source plane is larger then threshold
      % then stop.
   
      BetaMean      = mean(BetaAllIm);
      % residuals in the source plane
      Dist_S        = [BetaAllIm(:,1) - BetaMean(1); BetaAllIm(:,2) - BetaMean(2)];
   
      %--- Method to use in calculating source plane RMS ---
      switch CheckSourcePosScat
       case 'StD'
          %--- Check StD of source position ---
          RMS_S = sqrt(mean(Dist_S.^2));
      
       case 'MaxMean'
          %--- Check Max distance between sources position and mean source position ---
          RMS_S = sqrt(max( (BetaAllIm(:,1) - BetaMean(1)).^2 + (BetaAllIm(:,2) - BetaMean(2)).^2 ));
      
       case 'MaxMedian'
          %--- Check Max distance between sources position and median source position ---
          BetaMedian = median(BetaAllIm);
          RMS_S = sqrt(max( (BetaAllIm(:,1) - BetaMedian(1)).^2 + (BetaAllIm(:,2) - BetaMedian(2)).^2 ));
      
       case 'no'
          %--- do nothing ---
          RMS_S = 1e10;
       otherwise
          error('Unknown CheckSourcePosScatter Option');
      end
      
      if (RMS_S>SourcePosScatCrit),
         %--- Scatter in source position is too large ---
         % model is probably invalid
      
         %%%--- RMS = ??
   
         switch StopIf_S_ScatterLarge
          case 1
             Stop = 1;
          case 0
             Stop = 0;
          otherwise
             error('Unknown StopIf_S_ScatterLarge Option');
         end
   
      else
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
      
         %------------------------------------------------------
         %--- Check the scatter starting with a single image ---
         %------------------------------------------------------
         % i.e., don't solve for beta (psuedo image plane solution).
         switch CheckSingleBeta
          case 'yes'
             [RMS_PI,Dist_PI,ImPos] = find_images_regions(BetaStart,ModelPars,ModelType,ObsImPos,Dls_Ds(Is));

          case 'no'
             %--- do nothing ---
             RMS_PI              = 0;
             Dist_PI             = zeros(Nim.*2,1).*NaN;
          otherwise
             error('Unknown CheckSingleBeta Option');
         end
   
         switch StopIf_PI_ScatterLarge
          case 1
             Stop = 1;
          case 0
             Stop = 0;
          otherwise
             error('Unknown StopIf_PI_ScatterLarge Option');
         end
      
         if (RMS_PI>PseImagePosScatCrit),
            %--- Scatter is too large ---
            Dist_I = zeros(Nim.*2,1).*NaN;
            BestBeta = [NaN, NaN];
         else
            %---------------------------------------------------------------------------------------
            %--- Solve for sorce position (Beta) that is minimize the scatter in the image plane ---
            %---------------------------------------------------------------------------------------
            MinimizeMethod = 'Direct';
            [BestBeta,RMS_I,Dist_I,ImPos]=beta_minimize(BetaStart,ModelPars,ModelType,ObsImPos,Dls_Ds(Is),MinimizeMethod);
            switch StopIf_I_ScatterLarge
             case 1
                Stop = 1;
             case 0
                Stop = 0;
             otherwise
                error('Unknown StopIf_I_ScatterLarge Option');
            end
   
         end
      end
   
      AllBestBeta= [AllBestBeta; BestBeta];
      AllDist_S  = [AllDist_S;  Dist_S];
      AllDist_PI = [AllDist_PI; Dist_PI];
      AllDist_I  = [AllDist_I;  Dist_I];
      % all corresponding image position errors
      if (size(ObsImPos,2)>2),
         AllErr_I   = [AllErr_I;   ObsImPos(:,3); ObsImPos(:,4)];
      else
         AllErr_I   = NaN;
      end
   else
      %--- Stop criterion was activated ---
      % do nothing
   end
end


%--------------------------
%--- Calculate RMS/Chi2 ---
%--------------------------
RMS_S  = std(AllDist_S);
RMS_PI = std(AllDist_PI);
if (size(ObsImPos,2)>2),
   %--- Calculate chi2 ---
   RMS_I = sum( (AllDist_I./AllErr_I).^2 );
   Dof   = length(AllDist_I);
else
   %--- Calculate RMS ---
   RMS_I = std(AllDist_I);
   Dof   = length(AllDist_I);
end


