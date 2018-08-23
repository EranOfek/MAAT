function [Dist_SP,RMS_SP,BetaMean,Norm]=s1plane_rms(ModelPars,ModelType,ObsImPos);
%--------------------------------------------------------------------------
% s1plane_rms function                                               glens
% Description: Given a model parameters of a lens, and the
%                          images position corresponding to a single source,
%  calculate the deflections
%                          (and magnification) and the sources position,
%        and the best normalization.
%                          The program calculates the rms of the sources
%                          position in the source plane.
%               assuming Dls_Ds = 1.
% Input  : - Model parameters (see calc_alpha.m).
%          - Vector of Model Type (see calc_alpha.m).
%          - Images position [ThetaX, ThetaY, ErrorX, ErrorY]
%            It is recomended that the first image in each list
%            will be the faintest image corresponds to each source.
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

Nsource       = 1;
AllBetaMean   = zeros(Nsource,2);
AllDist_S     = zeros(0,1);

Nim      = size(ObsImPos,1);
%------------------------------------------------------------------
%--- Calculate deflection and magnification for images position ---
%------------------------------------------------------------------
[AlphaX,AlphaY] = calc_alpha(ObsImPos(:,1),ObsImPos(:,2),ModelPars,ModelType);
%--- Calculate source position corresponding to each image position ---

% fit source position and normalization using least squares
% [BetaX, BetaY, Norm]
H = [[ones(Nim,1); zeros(Nim,1)], [zeros(Nim,1); ones(Nim,1)], [AlphaX; AlphaY]];
ThetaVec = [ObsImPos(:,1); ObsImPos(:,2)];
Par = H\ThetaVec;
BetaMean  = [Par(1), Par(2)];
Norm      = Par(3);
%Resid     = H*Par - ThetaVec;
%RMS_Theta = std(Resid);   % equivalent to RMS

% Beta = Theta - Norm*Alpha
BetaX    = ObsImPos(:,1) - Norm.*AlphaX;
BetaY    = ObsImPos(:,2) - Norm.*AlphaY;
Dist_SP  = [BetaX - BetaMean(1); BetaY - BetaMean(2)];
RMS_SP   = std(Dist_SP);

