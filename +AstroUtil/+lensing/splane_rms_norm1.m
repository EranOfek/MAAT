function [RMS_S,Dist_S,AllBetaMean]=splane_rms_norm1(Norm,ModelPars,ModelType,ImagesCell,Dls_Ds);
%--------------------------------------------------------------------------
% splane_rms_norm1 function                                          glens
% Description: Given a model parameters of a lens, and the
%                          images position, calculate the deflections
%                          (and magnification) and the sources position.
%                          The program calculates the rms of the sources
%                          position in the source plane.
%                          ================================================
%                          This program is a duplicate of splane_rms.m
%                          The major difference is that this program get
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
% Output : - RMS in the source plane.
%          - Residual of each source corresoinding to an image, relative to the
%            mean position of the source.
%          - Mean source position for each set of images.
% Tested : Matlab 7.0
%     By : Eran Ofek          April 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See Also: splane_rms.m
%--------------------------------------------------------------------------
ColNorm     = 7;

CheckSourcePosScat     = 'StD';             % {'StD' | 'MaxMean'}

%--- noramlized Norm according to first line ---
NormalizedNorm        = ModelPars(:,ColNorm)./ModelPars(1,ColNorm);
%--- Set the Normalizations for all lines ---
ModelPars(:,ColNorm)  = Norm .* NormalizedNorm;



Nsource       = length(ImagesCell);
AllBetaMean   = zeros(Nsource,2);
AllDist_S     = zeros(0,1);
for Is=1:1:Nsource,
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
   
   AllDist_S           = [AllDist_S;  Dist_S];
   AllBetaMean(Is,:)   = BetaMean;
end


%-----------------------------------------------------
%--- Method to use in calculating source plane RMS ---
%-----------------------------------------------------
switch CheckSourcePosScat
 case 'StD'
    %--- Check StD of source position ---
    RMS_S = sqrt(mean(AllDist_S.^2));
      
 case 'MaxMean'
    %--- Check Max distance between sources position and mean source position ---
    RMS_S = max( AllDist_S );
      
 otherwise
    error('Unknown CheckSourcePosScatter Option');
end


