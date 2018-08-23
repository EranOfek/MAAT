function Res=fit_ellipse(X,Y,ParFlag)
% 2D ellipse fitting
% Package: AstroUtil.binary
% Description: Fit an ellipse of a form:
%              aX^2 + bY^2 + cX + dY +eX*Y = 1 to set of 2-D points.
% Input  : - Vector of X points.
%          - Vector of Y points.
%          - Vector of 5 logical flags indicating which ellipse parameters
%            to fit.
%            The logical flags corresponds to the [a, b, c, d, e]
%            parameters of the equation aX^2 + bY^2 + cX + dY +eX*Y = 1.
%            Default is true(1,5).
% Output : - Structure with the following fields:
%            .Par   - best fit parameters
%            .Resid - Vector of best fit residuals.
%            .Npar  - Number of free parameters.
%            .Ndof  - Number of degrees of freedom.
%            .RMS   - RMS.
%     By : Eran O. Ofek                    Dec 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=fit_ellipse(X,Y);
% Reliable: 2


if (nargin<3)
    ParFlag = true(1,5);
end

N = numel(X);

% The design matrix
H = [X(:).^2, Y(:).^2, X(:), Y(:), X(:).*Y(:)];
YY = ones(N,1);
H = H(:,ParFlag);
Res.Par = H\YY;

Res.Resid = H*Res.Par - 1;
Res.Npar  = numel(Res.Par);
Res.Ndof  = N - Res.Npar;
Res.RMS   = std(Res.Resid);

    





