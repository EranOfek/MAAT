function [Out,ParErr]=ls_conjgrad(H,Y,ErrY,Fun,varargin)
% Solve a linear least squares problem using the conjugate gradient method.
% Package: Util.fit
% Description: Solve a linear least squares problem using the conjugate
%              gradient method.
% Input  : - The design matrix.
%          - Vector of observations.
%          - Scalar of error, or vector of errors,
%            or a covariance matrix.
%          - One of the following functions to use:
%            'pcg'  - Preconditioned conjugate gradients method.
%            'cgs'  - Conjugate gradients squared method
%            'bicg' - Biconjugate gradients method
%            'bicgstab' - Biconjugate gradients stabilized method
%            'bicgstabl' - Biconjugate gradients stabilized(l) method
%            'qmr'  - Quasi-minimal residual method
%            'minres' - Minimum residual method
%            'gmres' - Generalized minimum residual method (with restarts)
%          * Arbitrary number of additional parameters to pass to
%            the minimizing function. Default is no parameters.
% Output : - A structure containing the following fields:
%            .Par    - Vector of free parameters.
%            .ParErr - Vector of free parameter errors.
%            .Resid  - Vector of residuals.
%            .Chi2   - \chi^2
%            .Npar   - Number of free parameters.
%            .Neq    - Number of equations.
%            .Dof    - Number of degrees of freedom.
%            .Cov    - Covariance matrix of parameters.
%            Or the vector of .Par (fitted parameters) of two output
%            arguments are requested.
%          - Errors in fitted parameters.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Nov 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=[0:1:100]'; Y=5+0.1.*X+randn(size(X)); H=[ones(size(X)), X];
%          Out=ls_conjgrad(H,Y,1,'pcg',1e-6,10)
% Reliable: 2
%--------------------------------------------------------------------------

%if (numel(ErrY)==1),
   % error is scalar
   A    = H.'*H;
   YY   = H.'*Y;
   Npar = size(H,1);
   %Cov  = inv(H'*inv(diag(ErrY.^2.*ones(Npar,1)))*H);
%end
% %else
%   if (min(size(ErrY))==1),
%       % error is vector
%       Cov = inv(H'*inv(diag(ErrY.^2))*H);
%       A   = Cov*H'*inv(diag(ErrY.^2));
%   else
%       % error is matrix (covariance matrix)
%       Cov = inv(H'*inv(ErrY)*H);
%       A   = Cov*H'*inv(ErrY);
%   end
% %end


Out.Par    = feval(Fun,A,YY,varargin{:});
Out.ParErr = nan(size(Out.Par));

if (nargout>1),
    ParErr = Out.ParErr;
    Out    = Out.Par;
end
    


%Out.ParErr = sqrt(diag(Cov));
% Out.Resid  = Y - H*Out.Par;
% Out.Chi2   = sum((Out.Resid./ErrY).^2);
% Out.Npar   = size(H,2);
% Out.Neq    = size(H,1);
% Out.Dof    = Out.Neq - Out.Npar;
%Out.Cov    = Cov;
