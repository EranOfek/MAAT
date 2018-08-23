function Res=fit_gauss1dp(X,Y,ErrY)
% SHORT DESCRIPTION HERE
% Package: Util.fit
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

if (nargin==2)
    ErrY = 1;
end

X0   = 4000;
LogY = log(Y);

Converge = false;
%while ~Converge
VecX0 = (4000:10:6000)';
Nx0   = numel(VecX0);
I = 0;
for I=1:1:Nx0; 
    X0 = VecX0(I);
    H = (X-X0).^2;
    Par = H\LogY;
    Sigma(I) = sqrt(0.5./Par);
    Resid = Y - exp(H*Par);
    Chi2(I) = sum(Resid.^2./Y);
    S(I) = std(Resid./LogY);
end
plot(VecX0,Chi2)
figure(2); plot(VecX0,Sigma)
Res = [];

