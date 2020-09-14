function [S2,Summary]=sysrem(Resid,Sigma,varargin)
% Apply the Tamuz et al. sysrem decomposition to a matrix of residuals
% Package: timeseries
% Description: Given a matrix of residuals (R_ij), of star i in image j,
%              iteratively decompose the matrix by minimizing
%              R'_ij = R_ij - C_i*A_j. The C and A vectors are free
%              parameters that are found by iterations.
%              Optionally, this can be applied iterativel several times -
%              e.g., R''_ij = R'_ij - C'_i * A'_j.
% Input  : - Matrix of residuals (R_ij), of star i in image j.
%          - Matrix of errors in residuals (Sigma_ij), of star i in image
%            j. If scalar then will be applied to all observations.
%            Default is 1.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'A' - Guess vector of A to be used in the first iteration.
%                  If scalar then will be duplicated into a vector.
%                  Default is 1.
%                  For Iter>1, 1 is used.
%            'C' - Guess vector of A to be used in the first iteration.
%                  If scalar then will be duplicated into a vector.
%                  Default is 0.
%                  For Iter>1, 0 is used.
%            'ReNormSigma' - Renormalize Sigma matrix (before first
%                  iterauion), such that the Chi^2 per dof will be 1.
%                  {true | false}. Default is true.
%            'ThreshDeltaS2' - Convergence threshold (in units of chi^2).
%                  Default is 1.
%            'Niter' - Number of iteration to apply sysrem. Default is 1.
% Output : - The resulted S2 value (\chi^2).
%          - A structure array of length (Niter+1), that contains the
%            results from each iteration, The last element is for the final
%            iteration.
%            Available fields:
%            'Resid' - Matrix of new residuals.
%            'A'     - Vector A
%            'B'     - Vector B
%            'S2'    - \chi^2
%            'rms'   - std of Resid matrix.
%            'Ndof   - Number of degrees of freedom.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [S2,Res]=timeseries.sysrem(R,1)
% Reliable: 2
%--------------------------------------------------------------------------

Def.Sigma = 1;
if (nargin<2)
    Sigma = Def.Sigma;
end

DefV.A                    = 1;
DefV.C                    = zeros;
DefV.ReNormSigma          = true;
DefV.ThreshDeltaS2        = 1;
DefV.Niter                = 1;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


[Nst,Nim] = size(Resid);
Nobs = Nst.*Nim;
Ndof = Nobs - Nst - Nim;

if (numel(InPar.A)==1)
    A = ones(1,Nim);
else
    A = InPar.A(:).';
end

if (numel(InPar.C)==1)
    C = ones(Nst,1);
else
    C = InPar.C(:);
end

S2_Orig = sum( (Resid./Sigma).^2,[1 2]);

if (InPar.ReNormSigma)
    Sigma = Sigma.*sqrt(S2_Orig./Ndof);

    S2_Prev = sum( (Resid./Sigma).^2,[1 2]);
else
    S2_Prev = S2_Orig;
end

K = 1;
Summary(K).Resid = Resid;
Summary(K).A     = A;
Summary(K).C     = C;
Summary(K).S2    = S2_Prev;
Summary(K).rms   = std(Resid(:));
Summary(K).Ndof  = Nobs;
 

for IterR=1:1:InPar.Niter
    K = K + 1;
    
    if (IterR>1)
        C = zeros(Nst,1);
        A = ones(1,Nim);
    end
    
    DeltaS2 = Inf;

    Iter = 0;
    while abs(DeltaS2)>InPar.ThreshDeltaS2
        Iter = Iter + 1;
        %InPar.Niter

        C = sum(Resid.*A ./ (Sigma.^2), 2) ./ sum( A.^2./Sigma.^2, 2);
        A = sum(Resid.*C ./ (Sigma.^2), 1) ./ sum( C.^2./Sigma.^2, 1);

        S2 = sum( ((Resid - C.*A)./Sigma).^2,[1 2]);
        DeltaS2 = S2 - S2_Prev;
        S2_Prev = S2;
    end

    Resid = Resid - C.*A;
    
    Summary(K).Resid = Resid;
    Summary(K).A     = A;
    Summary(K).C     = C;
    Summary(K).S2    = S2;
    Summary(K).rms   = std(Resid(:));
    Summary(K).Ndof  = Summary(K-1).Ndof - Nst - Nim;
    
end