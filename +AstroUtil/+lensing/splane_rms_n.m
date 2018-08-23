function [AllDist_S,RMS_S,AllBetaMean,BestNorm]=splane_rms_n(ModelPars,ModelType,ImagesCell,Dls_Ds);
%--------------------------------------------------------------------------
% splane_rms_n function                                              glens
% Description:       Given a model parameters of a lens, and the
%                          images position, calculate the deflections
%                          (and magnification) and the sources position,
%        and the best normalization.
%                          The program calculates the rms of the sources
%                          position in the source plane.
% Input  : - Model parameters (see calc_alpha.m).
%          - Vector of Model Type (see calc_alpha.m).
%          - Images position cell array:
%            Each cell contains [ThetaX, ThetaY, ErrorX, ErrorY]
%            of the images of a single source [pixels units!].
%            It is recomended that the first image in each list
%            will be the faintest image corresponds to each source.
%          - Vector of Dls/Ds for each source.
% Output : - Residual of each source corresoinding to an image, relative to the
%            mean position of the source.
%          - RMS in the source plane.
%          - Mean source position for each set of images.
%          - Best fit normalization (multiply the input normalization),
%            for each set of images.
% Tested : Matlab 7.0
%     By : Eran Ofek          April 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See Also: splane_rms_norm.m
%--------------------------------------------------------------------------

CheckSourcePosScat     = 'StD';             % {'StD' | 'MaxMean'}

Nsource       = length(ImagesCell);
AllBetaMean   = zeros(Nsource,2);
AllDist_S     = zeros(0,1);

AlphaX = cell(Nsource);
AlphaY = cell(Nsource);
for Is=1:1:Nsource,
   ObsImPos = ImagesCell{Is};
   Nim      = size(ObsImPos,1);
   %------------------------------------------------------------------
   %--- Calculate deflection and magnification for images position ---
   %------------------------------------------------------------------
   %[AlphaX,AlphaY,A_11,A_22,A_12] = calc_alpha(ObsImPos(:,1),ObsImPos(:,2),ModelPars,ModelType);
   [AlphaX{Is},AlphaY{Is}] = calc_alpha(ObsImPos(:,1),ObsImPos(:,2),ModelPars,ModelType);
   %--- Calculate source position corresponding to each image position ---


   % find the best normalization and source position using least squares (no weights).
   % [BetaX, BetaY, Norm]
   H = [[ones(Nim,1); zeros(Nim,1)], [zeros(Nim,1); ones(Nim,1)], Dls_Ds(Is).*[AlphaX; AlphaY]];
   ThetaVec = [ObsImPos(:,1); ObsImPos(:,2)];
   Par = H\ThetaVec;
   BetaMeanX = Par(1);
   BetaMeanY = Par(2);
   Norm      = Par(3);

   Resid = H*Par - ThetaVec;

   % Beta = Theta - Norm*Alpha
   BetaX   = ObsImPos(:,1) - Norm.*AlphaX;
   BetaY   = ObsImPos(:,2) - Norm.*AlphaY;
   Dist_SX = BetaX - BetaMeanX;
   Dist_SY = BetaY - BetaMeanY;

   AllDist_S           = [AllDist_S; Dist_SX; Dist_SY];
   AllBetaMean(Is,:)   = [BetaMeanX, BetaMeanY];
   BestNorm(Is)        = Norm;

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




