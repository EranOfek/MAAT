function ModelPars=lensmodel_solvenorm(ModelPars,ModelType,ImagesCell,Dls_Ds,Solve_IP);
%---------------------------------------------------------------------------
% lensmodel_solvenorm function                                        glens
% Description:    Given a lens model and images position
%                                solve for the best fit normalization
%                                (7th column parameter in ModelPars).
%                                The program (can) first find the best
%                                fit solution in the source plane and then
%                                use this to do the fit in the image plane.
% Input  : - Model parameters (see calc_alpha.m).
%          - Vector of Model Type (see calc_alpha.m).
%          - Images position cell array:
%            Each cell contains [ThetaX, ThetaY, ErrorX, ErrorY]
%            of the images of a single source [pixels units!].
%            It is recomended that the first image in each list
%            will be the faintest image corresponds to each source.
%          - Vector of Dls/Ds for each source.
%          - Solve for normalization using image plane minimization {1 | 0},
%            default is 1 (otherwise use source plane minimization).
% Output : - Model parameters (see calc_alpha.m) in which the normalization
%            (7th column) is the best fit source plane solution.
% Tested : Matlab 7.0
%     By : Eran O. Ofek            April 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------
ColNorm = 7;                   % Normalization column in ModelPars
CalcAlphaOnce = 1;             % calc deflection field once

Solve_SP = 1;                  % {1 | 0}
Method_SP = 'fminsearch';      % {'fminsearch'}
if (nargin==4),
   Solve_IP = 1;                  % {1 | 0}
end
Method_IP = 'fminsearch';      % {'fminsearch'}

switch Solve_SP
 case 1
    %---------------------------------------------------
    %--- Find best fit normalization in source plane ---
    %---------------------------------------------------

    Norm0           = ModelPars(1,ColNorm);   % first Guess for Normalization
    %--- noramlized Norm according to first line --
    NormalizedNorm  = ModelPars(:,ColNorm)./Norm0;
    Options         = optimset('fminsearch');


    switch Method_SP
     case 'fminsearch'
        %-----------------------------
        %--- Use MATLAB fminsearch ---
        %-----------------------------

        switch CalcAlphaOnce
         case 0
            Norm  = fminsearch('splane_rms_norm1',Norm0,Options,ModelPars,ModelType,ImagesCell,Dls_Ds);
            %--- Final source plane best fit normalization ---
            ModelPars(:,ColNorm) = Norm .* NormalizedNorm;

         case 1
            Ns         = length(ImagesCell);
            CellAlpha  = cell(Ns,1);
            for Is=1:1:Ns,
               [AlphaX,AlphaY] = calc_alpha(ImagesCell{Is}(:,1),ImagesCell{Is}(:,2),ModelPars,ModelType);
               CellAlpha{Is} = [AlphaX, AlphaY];
            end
            Norm0 = 1;
            Norm  = fminsearch('splane_rms_norm',Norm0,Options,CellAlpha,ImagesCell,Dls_Ds);
            ModelPars(:,ColNorm) = ModelPars(:,ColNorm).*Norm;
         otherwise
            error('Unknown CalcAlphaOnce Option');
        end

     otherwise
        error('Unknonw Method_SP Option');
    end


 case 0
    % do nothing
 otherwise
    error('Unknown Solve_SP Options');
end

switch Solve_IP
 case 1
    %---------------------------------------------------
    %--- Find best fit normalization in source plane ---
    %---------------------------------------------------

    Norm0           = ModelPars(1,ColNorm);   % first Guess for Normalization
    %--- noramlized Norm according to first line --
    NormalizedNorm  = ModelPars(:,ColNorm)./Norm0;

    switch Method_IP
     case 'fminsearch'
        %-----------------------------
        %--- Use MATLAB fminsearch ---
        %-----------------------------
        Norm  = fminsearch('iplane_rms_norm',Norm0,Options,ModelPars,ModelType,ImagesCell,Dls_Ds);

     otherwise
        error('Unknonw Method_SP Option');
    end

    %--- Final source plane best fit normalization ---
    ModelPars(:,ColNorm) = Norm .* NormalizedNorm;

 case 0
    % do nothing
 otherwise
    error('Unknown Solve_SP Options');
end
