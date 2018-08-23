function [BestBeta,RMS,Dist,ImPos]=beta_minimize(BetaStart,ModelPars,ModelType,ObsImagesPos,DlsDs,MinimizeMethod);
%-----------------------------------------------------------------------------
% beta_minimize function                                                glens
% Description: Given a model and approximate source position,
%                          find the best source position that minimize the
%                          residuals in the image plane.
% Input  : - Source position to start with [BetaX, BetaY].
%          - ModelPars (see calc_alpha.m).
%          - ModelType (see calc_alpha.m).
%          - Observed images position [ThetaX, ThetaY],
%            or [X, Y, ErrX, ErrY] (in that case return chi2 instead of RMS).
%          - Dls/Ds
%          - Minimization method:
%            'Direct'   - direct search using fminsearch.
% Output : - Best source position that minimize the residuals in the image
%            plane [BetaX, BetaY].
%          - RMS or chi2.
%          - Difference between observed and calculated (best) image position,
%            x followed by y.
%          - Calculated (best) images position.
% Tested : Matlab 7.0
%     By : Eran O. Ofek              April 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------------------

import AstroUtil.lensing.*

switch MinimizeMethod
 case 'Direct'
    %---------------------------------------------------
    %--- Direct minimization using MATLAB fminsearch ---
    %---------------------------------------------------
    Options  = optimset('fminsearch');
    Options  = optimset(Options,'TolFun',1e-3,'TolX',1e-3,'MaxFunEvals',10,'MaxIter',10);


    BestBeta = fminsearch('find_images_regions',BetaStart,Options,ModelPars,ModelType,ObsImagesPos,DlsDs);
    %--- RMS for BestBeta ---
    [RMS,Dist,ImPos] = find_images_regions(BestBeta,ModelPars,ModelType,ObsImagesPos,DlsDs);

 otherwise
    error('Unknown MinimizeMethod Option');
end


