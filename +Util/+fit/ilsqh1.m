function [A,B]=ilsqh1(X,Sigma,Norm)
%--------------------------------------------------------------------------
% ilsqh1 function                                                   FitFun
% Description: Solve the homogenized equation b_{j}*X_{ij}-a_{i}=0 in the
%              least square sense, with the constraint sum(a_{i})=N.
%              Where X_{ij} are the observables and
%              a_{i} and b_{j} are the unknowns.
% Input  : - Matrix X_{ij}.
%          - Matrix of errors \sigma_{ij}.
%          - Constraint for sum(a_{i}). Default is 1.
% Output : - Vector of a_{i}.
%          - Vector of b_{j}.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Oct 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: I=(-10:1:10); X=exp(-I.^2); X=repmat(X,10,1);
%          X = bsxfun(@times,X,(1:1:10)')+randn(size(X)).*0.01;
%          [A,B]=ilsqh1(X,ones(size(X)),1)
% Reliable: 
%--------------------------------------------------------------------------

if (nargin<3),
    Norm = 1;
end

[Na,Nb] = size(X);

% first guess for a and b
A = ones(Na,1)./Na.*Norm;
B = ones(Nb,1);

SumXS = sum((X./Sigma).^2,1);
SumS  = sum((1./Sigma).^2,2);

for Iiter=1:1:20,
    B = sum(bsxfun(@times,A,X)./(Sigma.^2),1)./SumXS;
    A = sum(bsxfun(@times,B,X)./Sigma.^2,2)./SumS;
    A = A./sum(A).*Norm;
    A
end




    