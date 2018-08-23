function [H,CalibH]=design_matrix_calib(P,Q,Sparse,Class)
%--------------------------------------------------------------------------
% design_matrix_calib function                                      FitFun
% Description: Build a design matrix for a linear system of the form:
%              A_ij = Z_i + M_j, where Z_i is of the length P, and M_j is
%              of length Q.
%              This is the design matrix for the basic calibration problem
%              (e.g., photometric calibration).
%              Note that the corresponding vector of observables contains
%              a concatenation the vectors of sources of each image.
% Input  : - P (e.g., number of images in the phot. calib. problem).
%            If method is 'index' and output is sparse then P and Q may be
%            vectors by which to populate the design matrix.
%          - Q (e.g., number of stars in the phot. calib. problem).
%          - A logical flag indicating if to return a sparse matrix.
%            Default is true.
%          - Output class (relevant only for non-sparse matrices).
%            Default is 'double'.
% Output : - The design matrix.
%          - The calibration block of the design matrix.
%            The rank of the design matrix (first output argument) is
%            P+Q-1. Hence, there is an additive freedom in the global ZP of
%            the solution. This block can concatenated to the design matrix
%            in order to force the global ZP.
%            (See Equation A4 in Ofek et al. 2011).
% Reference: Appendix A in Ofek et al. 2011 (ApJ, 740, 65).
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: H=Util.fit.design_matrix_calib(300,400);
% Reliable: 2
%--------------------------------------------------------------------------

Def.Sparse   = true;
Def.Class    = 'double';
if (nargin==2)
    Sparse   = Def.Sparse;
    Class    = Def.Class;
elseif (nargin==3)
    Class    = Def.Class;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments - Input (Q.M,[Sparse,Class])');
end

if (Sparse)
    Method = 'index';
    switch lower(Method)
        case 'index'
            % see comment about construction - below:
%             PisNum = false;
%             QisNum = false;
%             if (numel(P)>1)
%                 PV = P;
%                 P  = numel(P);
%                 PisNum = true;
%             end
%             if (numel(Q)>1)
%                 QV = Q;
%                 Q  = numel(Q);
%                 QisNum = true;
%             end
            AllInd1 = bsxfun(@plus,(1:1:Q) , (0:(P.*Q+Q):(P.*Q.*P)).' );
            AllInd2 = bsxfun(@plus,(0:Q:P.*Q-1), (1:(P.*Q+1):P.*Q.*Q)');
            AllInd2 = AllInd2 + P.*P.*Q;
            [I,J]=ind2sub([Q.*P, Q+P],[AllInd1(:);AllInd2(:)]);
            H = sparse(I,J,1);
%             % populate with specific numbers
%             if (PisNum)
%                 H(AllInd1(:)) = repmat(PV,Q,1).';
%             end
%             if (QisNum)
%                 H(AllInd2(:)) = repmat(QV,P,1).';
%             end
        case 'block'
            H = sparse([],[],[],Q.*P,P+Q,2.*P.*Q);
            Diag = diag(ones(Q,1));
            for Ip=1:1:P
                Ir = (Ip-1).*P + (1:1:Q);
                H(Ir,Ip) = 1;
                H(Ir,P+1:P+Q) = Diag; %diag(ones(Q,1));
            end
        otherwise
            error('Unknown Method option');
    end
                    
    
else
    H = zeros(P.*Q,P+Q,Class);
    % first block matrices
    % Columns 1 to P
    % Each sub matrix has Q rows, and in total there are P sub matrices
    % P*Q rows
    AllInd1 = bsxfun(@plus,(1:1:Q) , (0:(P.*Q+Q):(P.*Q.*P)).' );
    H(AllInd1) = 1;

    % second block matrices:
    % Columns P+1 to P+Q
    % These are a set of QxQ identity matrices one below the other...
    % Note that the output is [J,I] and NOT [I,J]:
    % (1:Q+1:Q.*Q) produce a series of Q identity matrices, each one of size QxQ
    AllInd2 = bsxfun(@plus,(0:Q:P.*Q-1), (1:(P.*Q+1):P.*Q.*Q)');
    AllInd2 = AllInd2 + P.*P.*Q;
    H(AllInd2) = 1;
end

if (nargout>1)
    % Set up the calibration block
    if (Sparse)
        K = P.*Q + (1:Q+1:Q.*Q);
        [I,J]=ind2sub([Q, Q+P],K);
        CalibH = sparse(I,J,1);
    else
        CalibH = [zeros(Q,P), diag(ones(Q,1))];
    end
end
        