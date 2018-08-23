function [BPar,BParErr,WPar,N1sig]=fitlin_ransac(X,Y,Sigma,Nsig,Nsim)
%--------------------------------------------------------------------------
% fitlin_ransac function                                            FitFun
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jan 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X = (1:1:1e6)';
%          Y= 2+ 0.95.*X + randn(size(X)).*sqrt(X);
%          Y(randi(1e6,8e5,1)) = rand(8e5,1).*1e6;
%          [Par]=fitlin_ransac(X,Y,sqrt(X));
% Reliable: 
%--------------------------------------------------------------------------

Def.Sigma = [];
Def.Nsig  = 1;
Def.Nsim  = 100;
if (nargin==2)
    Sigma = Def.Sigma;
    Nsig  = Def.Nsig;
    Nsim  = Def.Nsim;
elseif (nargin==3)
    Nsig  = Def.Nsig;
    Nsim  = Def.Nsim;
elseif (nargin==4)
    Nsim  = Def.Nsim;
elseif (nargin==5)
    % do nothing
else
    error('Illegal number of imput arguments');
end

N = numel(X);

%RandI = randi(N,Nsim,2);

H     = [ones(2,1), zeros(2,1)];
Par   = zeros(Nsim,2);
Dx    = zeros(Nsim,1);
N1sig = zeros(Nsim,1);

for Isim=1:1:Nsim
    RandI = randperm(N,2);
    %H(:,2)      = X(RandI(Isim,:));
    H(:,2)      = X(RandI);
    %Par(Isim,:) = (H\Y(RandI(Isim,:))).';
    %Dx(Isim)    = abs(diff(X(RandI(Isim,:))));
    Par(Isim,:) = (H\Y(RandI)).';
    Dx(Isim)    = abs(diff(X(RandI)));
    
    N1sig(Isim) = sum(abs(Y - ( Par(Isim,1) + Par(Isim,2).*X))<Nsig.*Sigma);
    % N1sig(Isim) = length(find(abs(Y - ( Par(Isim,1) + Par(Isim,2).*X))<Sigma));  % slow
end

%N1sig = sum(bsxfun(@lt,abs(bsxfun(@minus,Y,bsxfun(@plus,bsxfun(@times,X,Par(:,2).'),Par(:,1).'))),Nsig.*Sigma)).';

%hist(N1sig,100)
[MaxN,MaxInd] = max(N1sig);
%Par = Par(MaxInd,:);
WPar(1) = Util.stat.wmedian(Par(:,1),1./N1sig);
WPar(2) = Util.stat.wmedian(Par(:,2),1./N1sig);

% best fit using best points
Flag  = abs(Y - ( WPar(1) + WPar(2).*X))<Nsig.*Sigma;
N     = sum(Flag);
X     = X(Flag);
Y     = Y(Flag);
Sigma = Sigma(Flag);
H     = [ones(N,1), X];
[BPar,BParErr] = lscov(H,Y,1./(Sigma.^2));
%std(H*BPar - Y)





