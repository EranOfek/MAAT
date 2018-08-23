function [Par,Par_Err,Cov,Chi2,Freedom,Resid]=fitpoly(X,Y,DelY,Deg,SigClip,PlotPar)
%0-------------------------------------------------------------------------
% fitpoly function                                                  FitFun
% Description: Linear least squares polynomial fitting. Fit polynomial of th
%              form: Y= a_0 + a_1*X + a_2*X^2 +...+ a_n*X^n
%              to set of data points. Return the parameters, thir errors,
%              the \chi^2 squars, and the covariance matrix.
%              This function is replaced by fitgenpoly.m
% Input  : - Column vector of the independent variable.
%          - Column Vector of the dependent variable.
%          - Vector of the std error in the dependent variable.
%            If only one value is given, the points
%            are taken to be with equal weight. and Std error
%            equal to the value given.
%            If two columns are given then the second column is taken
%            as the error in the independent variable,
%            and the problem is solve iteratively, starting with b=0.
%            (Num. Rec. chapter 15).
%          - Degree of polynomial. (Default is 1).
%          - Sigma-Clipping (default is NaN, for no clipping).
%          - Vector of plot's control characters.
%            If argument is given then X vs. Y graph is plotted.
%            If equal to empty string (e.g. '') then plot X vs. Y
%            with red fitted function line and yellow circs for
%            the observations.
%            If one or two character are given then the first character
%            is for the observations sign and the second for the fitted
%            function line.
%            If third character is given then histogram of resdiual
%            is plotted. when the third character should contain the
%            number of bins.
% Output : - Fitted parameters [a0,a1,...]
%          - Fitted errors in the parameters [Da0,Da1,...]
%          - The covariance matrix.
%          - Chi2 of the fit.
%          - Degrees of freedom.
%          - The Y axis residuals vector.
% See also: fitslope.m, fit_lin.m
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Mar 1995
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 1
%--------------------------------------------------------------------------

MaxNIter = 5;   % maximum number of sigma-clipping iterations
if (nargin<4),
   Deg = 1;
   SigClip = NaN;
elseif (nargin<5),
   SigClip = NaN;
else
   % do nothing
end

N_X  = length(X);
N_Y  = length(Y);
N_DY = size(DelY,1);
N_DX = size(DelY,2);
if (N_DX==2),
   DelX = DelY(:,2);
   DelY = DelY(:,1);
else
   DelX = zeros(N_DY,1);
end
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
      X    = X(SCInd);
      Y    = Y(SCInd);
      DelY = DelY(SCInd);
      DelX = DelX(SCInd);
      
      N_X  = length(X);
      N_Y  = length(Y);
      N_DY = length(DelY);  
   end

   % degree of freedom
   Freedom = N_X - (Deg + 1);
   
   % building the H matrix
   H = zeros(N_X,Deg+1);
   H(:,1) = ones(N_X,1);
   for Ind=2:1:Deg+1,
      H(:,Ind) = X.^(Ind-1);
   end
   
   % building the Covariance matrix
   B       = 0;
   DelXY   = sqrt(DelY.^2 + (B.*DelX).^2);
   
   % Old - Memory consuming
   V       = diag(DelXY.^2);
   Cov     = inv(H'*inv(V)*H);
   Par     = Cov*H'*inv(V)*Y;
   Par_Err = sqrt(diag(Cov));

%%   [Par,Par_Err] = lscov(H,Y,DelXY.^2)   
   
   %'Number of degree of freedom :', Freedom
   Resid = Y - H*Par;
   Chi2  = sum((Resid./DelXY).^2);

   %Chi2/Freedom
   %sqrt(2/Freedom)
   
   NB_Iter = 1;
   while (abs(B-Par(2))>Par_Err(2).*1e-8),
      % iterate with DelX (B~=0)
      B = Par(2);   
      % building the Covariance matrix
      DelXY   = sqrt(DelY.^2 + (B.*DelX).^2);
      V       = diag(DelXY.^2);
      % Old - Memory consuming
      Cov     = inv(H'*inv(V)*H);
      Par     = Cov*H'*inv(V)*Y;
      Par_Err = sqrt(diag(Cov));
      %'Number of degree of freedom :', Freedom
      Resid = Y - H*Par;
      Chi2  = sum((Resid./DelXY).^2);
      NB_Iter = NB_Iter + 1;
   end   
   
end

if (nargin==6),
   % plot results
   length(PlotPar);
   if (length(PlotPar)==0),
      PlotPar(1) = 'o';
      PlotPar(2) = 'r';
   end
   if (length(PlotPar)==1),
      PlotPar(2) = 'r';
   end
   figure(1);
   plot(X,Y,PlotPar(1));
   hold on;
   plot(X,H*Par,PlotPar(2));
   xlabel('X');
   ylabel('Y');
   hold off;
   if (length(PlotPar)==3),
      % plot histogram of residuals
      figure(2);
      [Hist_X,Hist_N]=realhist(sort(abs(Resid)),str2num(PlotPar(3)),[0,max(abs(Resid)).*1.0001]);
      bar(Hist_X,Hist_N);
      axis([0,max(abs(Resid)).*1.0001,0,max(Hist_N)+1]);
      xlabel('X');
      ylabel('Number');
   end
end


%errorxy([X,Y,DelY],[1 2 3],'.');
%hold on;
%plot(X,H*Par,'r');


fprintf(1,'\n Number of iterations : %d \n',Iter);


