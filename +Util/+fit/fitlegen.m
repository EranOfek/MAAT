function [Par,Par_Err,Cov,Chi2,Freedom,Resid]=fitlegen(X,Y,DelY,Deg,PlotPar);
%------------------------------------------------------------------------------
% fitlegen function                                                     FitFun
% Description: Fit Legendre polynomials to data, where the fitted function has
%              the form: Y= a_0*L_0(X) + a_1*L_1(X) +...+ a_n*L_n(X)
%              This function is replaced by fitgenpoly.m
% Input  : - Column vector of the independent variable.
%          - Column Vector of the dependent variable.
%          - Vector of the std error in the dependent variable.
%            If only one value is given then, points
%            are taken to be with equal weight. and Std error
%            equal to the value given.
%          - Degree of Legendre polynomial. (Default is 1).
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
%          - \chi^2 of the fit.
%          - Degrees of freedom.
%          - The Y axis residuals vector.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                      June 1998
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

if (nargin<4),
   Deg = 1;
end
N_X  = length(X);
N_Y  = length(Y);
N_DY = length(DelY);
if (N_X~=N_Y),
   error('X and Y must have the same length');
end
if (N_X~=N_DY),
   if (N_DY==1)
      % take equal weights
      DelY = DelY.*ones(N_X,1);
   else
      error('Y and DelY must have the same length');
   end
end

% degree of freedom
Freedom = N_X - (Deg + 1);

% building the H matrix
AssocLegD = 1; % (associated Legendre of order 0)
H = zeros(N_X,Deg+1);
H(:,1) = legendre(0,X);
for Ind=2:1:Deg+1,
   TempLegen = legendre(Ind-1,X);
   H(:,Ind) = TempLegen(AssocLegD,:)';
end

% building the Covariance matrix
V = diag(DelY.^2);

% Old - Memory consuming
Cov     = inv(H'*inv(V)*H);
Par     = Cov*H'*inv(V)*Y;
Par_Err = sqrt(diag(Cov));


%'Number of degree of freedom :', Freedom
Resid = Y - H*Par;
Chi2 = sum((Resid./DelY).^2);

%Chi2/Freedom
%sqrt(2/Freedom)

%plot(X,Y)
%hold on;
%plot(X,H*Par,'o');

if (nargin==5),
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




