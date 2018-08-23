function [Par,ErrPar]=cosmo_pars(ParType)
% Cosmological parameters as measured by various experiments.
% Package: AstroUtil.cosmo
% Description: Return the cosmological parameters as measured by
%              various experiments.
% Input  : - Cosmological parameter source, options are:
%            'wmap3'    : WMAP+SNLS (Spergel et al. 2007).
%            'wmap5'    : WMAP+SN+BO (Komatsu et al. 2008).
%            'wmap9'    : WMAP9+BAO+H0 (Hinshaw et al. 2013).
%            'planck'   : planck full mission (Ade et al. 2015).
%            or alternatively a three or four element vector
%            containing: [H0, OmegaM, OmegaL OmegaRad], where
%            default for OmegaRad is 0.
%            If this parameter is a structure, return the structure
%            as output.
% Output : - Structure containing cosmological parameters.
%          - Structure containing cosmological parameters errors
%            The tructure is identical to the parameter structure,
%            but containing two values for each parameter
%            [1sigma lower error, 1sigma upper error].
% Tested : Matlab 7.0
%     By : Eran O. Ofek                     March 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [Par,ErrPar]=cosmo_pars('wmap3');
% Reliable: 2
%--------------------------------------------------------------------------

if (isstruct(ParType)==1),
   Par    = ParType;
   ErrPar = NaN;
else
   if (ischar(ParType)==0),
      % vector
      Par.H0        = ParType(1);
      Par.OmegaM    = ParType(2);
      Par.OmegaL    = ParType(3);
      if (length(ParType)==3),
         Par.OmegaRad = 0;
      else
         Par.OmegaRad = ParType(4);
      end
      ErrPar = NaN;
   else
      % string   
      switch lower(ParType)
       case 'wmap3'
          %------------------------------------------
          % WMAP + SNLS colsmological parameters
          % Spergel et al. (2007; astro-ph/0603449)
          %------------------------------------------
          
          % Hubble parameter
          Par.H0         = 70.4;
          ErrPar.H0      = [-1.6 +1.5];
          
          Par.OmegaM     = 0.268;
          ErrPar.OmegaM  = [-0.018 +0.018];
          
          Par.OmegaB     = 0.03105;
          ErrPar.OmegaB  = [NaN NaN];
          
          Par.OmegaK     = 0.986;
          ErrPar.OmegaK  = [-0.017 +0.017];
          
          Par.OmegaL     = 0.716;
          ErrPar.OmegaL  = [-0.055 +0.055];
      
          % Omega radiation
          Par.OmegaRad    = 0.0;
          ErrPar.OmegaRad = [-0.0 0.0];
          
          % Equation of state (quintessence)
          Par.W            = -1; %OmegaM_H2    = 0.1324;
          ErrPar.W         = [-0.0 0.0];

          % OmegaM x h^2
          Par.OmegaM_H2    = 0.1324;
          ErrPar.OmegaM_H2 = [-0.0041 +0.0042];
          
          % OmegaB x h^2
          Par.OmegaB_H2    = 0.02186;
          ErrPar.OmegaB_H2 = [-0.00068 +0.00068];
          
          % Optical depth
          Par.Tau        = 0.073;
          ErrPar.Tau     = [-0.028 +0.027];
          
          % Power-law slope
          Par.Ns         = 0.947;
          ErrPar.Ns      = [-0.015 +0.015];
          
          % sigma8
          Par.Sigma8     = 0.776;
          ErrPar.Sigma8  = [-0.032 +0.031];

       case 'wmap5'
          %--------------------------------------------
          % WMAP5 + SNLS + BO colsmological parameters
          % Komatsu et al. (2008; astro-ph/0803.0547)
          %--------------------------------------------
          
          % Hubble parameter
          Par.H0         = 70.1;
          ErrPar.H0      = [-1.3 +1.3];
          
          Par.OmegaM     =   0.268;
          ErrPar.OmegaM  = [-0.018 +0.018];
          
          Par.OmegaB     = 0.0462;
          ErrPar.OmegaB  = [-0.0015 +0.0015];
          
          Par.OmegaK     =  1;
          ErrPar.OmegaK  = [NaN NaN];
          
          Par.OmegaL     = 0.721;
          ErrPar.OmegaL  = [-0.015 +0.015];
      
          % Omega radiation
          Par.OmegaRad    = 0.0;
          ErrPar.OmegaRad = [-0.0 0.0];
          
          % Equation of state (quintessence)
          Par.W            = -1;
          ErrPar.W         = [-0.0 0.0];

          % OmegaM x h^2
          Par.OmegaM_H2    = 0.1369;
          ErrPar.OmegaM_H2 = [-0.0037 +0.0037];
          
          % OmegaB x h^2
          Par.OmegaB_H2    = 0.02265;
          ErrPar.OmegaB_H2 = [-0.00059 +0.00059];
          
          % Optical depth
          Par.Tau        = 0.084;
          ErrPar.Tau     = [-0.016 +0.016];
          
          % Power-law slope
          Par.Ns         = 0.960;
          ErrPar.Ns      = [-0.013 +0.014];
          
          % sigma8
          Par.Sigma8     = 0.817;
          ErrPar.Sigma8  = [-0.016 +0.016];

          % Age
          Par.Age        = 13.73;          % [Gyr]
          ErrPar.Age     = [-0.12 0.12];   % [Gyr]

       case 'wmap9'
          %--------------------------------------------
          % WMAP9 + BAO + H0 colsmological parameters
          % Hinshaw et al. (2013)
          %--------------------------------------------
          
          % Hubble parameter
          Par.H0         = 69.33;
          ErrPar.H0      = [-0.88 +0.88];
          Par.OmegaM     =   0.2408+0.0471;
          ErrPar.OmegaM  = [-0.009 +0.009];
          
          Par.OmegaB     = 0.0471;
          ErrPar.OmegaB  = [-0.001 +0.001];
          
          Par.OmegaK     =  1;
          ErrPar.OmegaK  = [NaN NaN];
          
          Par.OmegaL     = 1- 0.2408-0.0471;
          ErrPar.OmegaL  = [-0.009 +0.009];
      
          % Omega radiation
          Par.OmegaRad    = 0.0;
          ErrPar.OmegaRad = [-0.0 0.0];
          
          % Equation of state (quintessence)
          Par.W            = -1;
          ErrPar.W         = [-0.0 0.0];

          % OmegaM x h^2
          %Par.OmegaM_H2    = 0.1369;
          %ErrPar.OmegaM_H2 = [-0.0037 +0.0037];
          
          % OmegaB x h^2
          Par.OmegaB_H2    = 0.02266;
          ErrPar.OmegaB_H2 = [-0.00043 +0.00043];
          
          % Optical depth
          Par.Tau        = 0.088;
          ErrPar.Tau     = [-0.013 +0.013];
          
          % Power-law slope
          Par.Ns         = 0.971;
          ErrPar.Ns      = [-0.01 +0.01];
          
          % sigma8
          Par.Sigma8     = 0.817;
          ErrPar.Sigma8  = [-0.016 +0.016];

          % Age
          Par.Age        = 13.75;          % [Gyr]
          ErrPar.Age     = [-0.085 0.085];   % [Gyr]
          
       case 'planck'
           % Planck full mission cosmological parameters
           % https://arxiv.org/abs/1502.01589
           % TT, TE, EE+low P, +lensing+ext
           
           Par.H0 = 67.74;
           ErrPar.H0 = [0.46 0.46];
           
           Par.OmegaM = 0.3089;
           ErrPar.OmegaM = [0.0062 0.0062];
           
           Par.OmegaRad = 0;
           
           Par.OmegaL = 0.6911;
           ErrPar.OmegaL = [0.0062 0.0062];
           
           Par.OmegaB_h2  = 0.02230;
           ErrPar.OmegaB_h2 = [0.00014 0.00014];
           
           Par.sigma8 = 0.8159;
           ErrPar.sigma8 = [0.0086 0.0086];
           
           Par.z_re = 8.8;
           ErrPar.z_re = [1.1 1.2];
           
           Par.Age = 13.799;
           ErrPar.Age = [0.021 0.021];
           
       otherwise
          error('Unknown ParType option');
      end
   end
end
