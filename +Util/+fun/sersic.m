function [Sigma,Ftot]=sersic(Pars,R)
%--------------------------------------------------------------------------
% sersic function                                                  General
% Description: Calculate a Sersic profile (Sersic 1968).
% Input  : - Row vector of parameters. If matrix, then each column
%            represent a parameter and each row represent a different
%            model. Parameters are [Sigma_e, R_e, n].
%            Sigma_e is the surface brightness at R_e.
%            R_e is the effective radius (at which half the ligt resides).
%            n is the Sersic index. 4 for de Vaucouleurs, 1 for
%            exponential, 0.5 for Gaussian.
%          - Column vector of radial coordinates. Default is (0:0.1:10).';
% Output : - Surface brighntness vector/matrix per radial position
%            and model. Different columns represents different models.
%          - Vector of total (2d integrated) flux in each model.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: Peng et al. (2002)
%            http://iopscience.iop.org/1538-3881/124/1/266/fulltext/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

if (nargin==1),
    R = (0:0.1:10).';
end

% Sigma = Sigma_e .* exp(-Kappa.*((R./R_e).^(1./n) - 1))
Sigma = Pars(:,1).* exp(-Kappa.*((R./Pars(:,2)).^(1./Pars(:,3)) - 1));
    

Sigma_e = 1;
R_e     = 1;
n       = 4;
Kappa   = 7.67;


Ftot = 2.*pi.*R_e.^2.*Sigma_e.*exp(Kappa).*n.*Kappa.^(-2.*n).*gamma(2.*n);