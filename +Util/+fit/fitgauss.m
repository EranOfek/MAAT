function [NewPar,NewPar_Err,Cov,Chi2,Freedom,Resid]=fitgauss(X,Y,DelY,SigClip);
%------------------------------------------------------------------------------
% fitgauss function                                                     FitFun
% Description: Fit a Gaussian function to data, where the Gaussian has the
%              form: Y = A * exp(-0.5.*((X-X0)./s).^2).
% Input  : - Column vector of the independent variable.
%          - Column Vector of the dependent variable.
%          - Vector of the std error in the dependent variable.
%            If only one value is given, then the points
%            are taken to be with equal weight, and Std error
%            equal to the value given.
%          - Sigma-Clipping (default is NaN, for no clipping).
% Output : - Fitted parameters [A,X0,s]
%          - Fitted errors in the parameters [DA,DX0,Ds]
%          - The covariance matrix.
%          - \chi^2 of the fit.
%          - Degrees of freedom.
%          - The Y axis residuals vector.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                   October 1996
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

MaxNIter = 5;   % maximum number of sigma-clipping iterations
if (nargin<4),
   SigClip = NaN;
else
   % do nothing
end

Deg  = 3;
NewY = log(Y); 

N_X  = length(X);
N_Y  = length(Y);
N_DY = length(DelY);
if (N_X~=N_Y),
   error('X and Y must have the same length');
end
if (N_X~=N_DY),
   if (N_DY==1),
      % take equal weights
      if (DelY<=0),
         error('DelY must be positive');
      else
         DelY = DelY.*ones(N_X,1);
      end
   else
      error('Y and DelY must have the same length');
   end
end

Resid = zeros(size(DelY));
if (isnan(SigClip)),
   MaxNIter = 1;
end

Iter = 0;
while (Iter<MaxNIter & (max(abs(Resid)>DelY | Iter==0))),
   Iter = Iter + 1;

   % sigma clipping
   if (isnan(SigClip)),
      % do not sigma clip
   else
      SCInd = find((abs(Resid)./(SigClip.*DelY))<1);  % find non-outlayers
%size(SCInd)
      X    = X(SCInd);
      Y    = Y(SCInd);
      NewY = NewY(SCInd);
      DelY = DelY(SCInd);

      N_X  = length(X);
      N_Y  = length(Y);
      N_DY = length(DelY);  
   end

   % degree of freedom
   Freedom = N_X - (Deg + 1);
   
   % building the H matrix
   H = zeros(N_X,Deg);
   H(:,1) = ones(N_X,1);
   for Ind=2:1:Deg,
      H(:,Ind) = X.^(Ind-1);
   end
   
   % building the Covariance matrix
   V = diag(DelY.^2);
   % Old - Memory consuming
   Cov     = inv(H'*inv(V)*H);
   Par     = Cov*H'*inv(V)*NewY;
   Par_Err = sqrt(diag(Cov));

   NewPar        = zeros(3,1);
   NewPar_Err    = zeros(3,1);
   NewPar(3)     = sqrt(-1./(2.*Par(3)));
   NewPar(2)     = NewPar(3).^2.*Par(2);
   NewPar(1)     = exp(Par(1) + 0.5.*NewPar(2).^2./(NewPar(3).^2));
   NewPar_Err(3) = (0.5./sqrt(2)).*Par_Err(3).*abs(Par(3)).^(-1.5);
   NewPar_Err(2) = sqrt((2.*NewPar(3).*Par(2).*NewPar_Err(3)).^2 + ...,
                        (NewPar(3).^2.*Par_Err(2)).^2);
   NewPar_Err(1) = sqrt((NewPar(1).*Par_Err(1)).^2 + ...,
                        (NewPar(1).*NewPar(2).^2.*NewPar(3).^(-3).*NewPar_Err(3)).^2 + ...,
                        (NewPar(1).*NewPar(2).*NewPar_Err(2)./(NewPar(3).^2)).^2);

                         
   
   %'Number of degree of freedom :', Freedom
   Resid = NewY - H*Par;
   Chi2  = sum((Resid./DelY).^2);

   %Chi2/Freedom
   %sqrt(2/Freedom)
end




%errorxy([X,Y,DelY],[1 2 3],'.');
%hold on;
%plot(X,H*Par,'r');


fprintf(1,'\n Number of iterations : %d \n',Iter);


