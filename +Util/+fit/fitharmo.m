function [Par,Par_Err,Cov,Chi2,Freedom,Par1,Resid]=fitharmo(X,Y,DelY,Har,Deg,PlotPar,File);
%------------------------------------------------------------------------------
% fitharmo function                                                     FitFun
% Description: Fit trignometric functions with their harmonies to data.
%              The fitted function is of the form:
%                      Y= a_1*sin(w1*t)     + b_1*cos(w1*t)   +
%                         a_2*sin(2*w1*t)   + b_2*cos(2*w1*t) + ...
%                         a_n*sin(n_1*w1*t) + b_n*cos(n_1*w1*t) + ...
%                         c_1*sin(w2*t)     + d_1*cos(w2*t) + ...
%                         s_0 + s_1*t + ... + s_n.*t.^n_s
%              Note that w is angular frequncy, w=2*pi*f, while the use
%              should indicate the frequency.
% Input  : - Column vector of the independent variable.
%          - Column Vector of the dependent variable.
%          - Vector of the std error in the dependent variable.
%            If only one value is given, the points
%            are taken to be with equal weight. and Std error
%            equal to the value given.
%          - matrix of harmonies to fit.
%            N*2 matrix, where N is the number of different frequncies.
%            Each row should contain two numbers, the first is the
%            frequency to fit and the second is the number of harmonies
%            of that frequncy to fit. If there is more then one row
%            then all the frequncies and their harmonics will be fitted
%            simoltanusly.
%          - Degree of polynomials to fit. (Default is 0).
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
%          - File name in which summary table of the fit will be written.
%            The summary table includes all information regarding the
%            fit parameters and Chi2 test.
% Output : - Fitted parameters [a_1,b_1,...,a_n,b_n,c_1,d_1,...,s_0,...]
%            The order of the parameters is like the order of the
%            freqencies matrix, and then the constant + linear terms.
%          - Fitted errors in the parameters [Da_1,Db_1,...]
%          - The covariance matrix.
%          - \chi^2 of the fit.
%          - Degrees of freedom.
%          - sine/cosine parameters in form od Amp. and phase (in fraction),
%            pairs of lines for [Amp, Amp_Err; Phase, Phase_Err]...
%            phase are given in the range [-0.5,0.5].
%          - The Y axis residuals vector.
% See also : fitharmonw.m
% Tested : Matlab 5.1
%     By : Eran O. Ofek                       May 1994
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

if (nargin<5),
   Deg = 1;
end
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

% number of parameters
N_Pars = Deg+1+2.*sum(Har(:,2));

% degree of freedom
Freedom = N_X - N_Pars;

% the size of the harmonies matrix
[Srow_Har,Scol_Har] = size(Har);
if (Scol_Har~=2),
   error('Number of columns in the harmonic freq. should be two');
end

% building the H matrix
H = zeros(N_X,N_Pars);
Counter = 0;
for I=1:1:Srow_Har,
   % run over number of harmonic per frequncy
   for J=1:1:Har(I,2),
      Counter = Counter + 1;
      H(:,Counter) = sin(2.*pi.*Har(I,1).*J.*X);
      Counter = Counter + 1;
      H(:,Counter) = cos(2.*pi.*Har(I,1).*J.*X);
   end
end
% add the constant term
Counter = Counter + 1;
H(:,Counter) = ones(N_X,1);
% add the linear terms
for I=1:1:Deg,
   Counter = Counter + 1;
   H(:,Counter) = X.^I;
end

% building the Covariance matrix
V = diag(DelY.^2);

% Old - Memory consuming
Cov     = inv(H'*inv(V)*H);
Par     = Cov*H'*inv(V)*Y;
Par_Err = sqrt(diag(Cov));

%'Number of degree of freedom :', Freedom
Resid = Y - H*Par;
Chi2  = sum((Resid./DelY).^2);

%Chi2/Freedom
%sqrt(2/Freedom)

if (nargin>5),
   'Plot Data'
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

   % ------------------------------------ New
   NewX = [min(X):((max(X)-min(X))./500):max(X)]';
   % building the newH matrix
   NewH = zeros(length(NewX),N_Pars);
   Counter = 0;
   for I=1:1:Srow_Har,
      % run over number of harmonic per frequncy
      for J=1:1:Har(I,2),
         Counter = Counter + 1;
         NewH(:,Counter) = sin(2.*pi.*Har(I,1).*J.*NewX);
         Counter = Counter + 1;
         NewH(:,Counter) = cos(2.*pi.*Har(I,1).*J.*NewX);
      end
   end
   % add the constant term
   Counter = Counter + 1;
   NewH(:,Counter) = ones(length(NewX),1);
   % add the linear terms
   for I=1:1:Deg,
      Counter = Counter + 1;
      NewH(:,Counter) = NewX.^I;
   end
   %----------
   plot(NewX,NewH*Par,PlotPar(2));
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

% write summary file
if (nargin==7),
   Fid = fopen(File,'w');
   fprintf(Fid,'         %s \n','fitharmo.m - summary file');
   fprintf(Fid,' %s \n',['Created by Eran O. Ofek   at : ',date]);
   fprintf(Fid,'\n');
   fprintf(Fid,'%s \n','F# H#     Frequncy    1/Frequncy        A       dA           B       dB           C       dC          Fi       dFi');

   F_Coun = 0;
   for I=1:1:Srow_Har,
      % run over number of harmonic per frequncy
      for J=1:1:Har(I,2),
         F_Coun = F_Coun + 1;
         C  = sqrt(Par(F_Coun).^2 + Par(F_Coun+1).^2);
         Fi = atan2(Par(F_Coun+1),Par(F_Coun));
         C_Err  = sqrt((Par(F_Coun).*Par_Err(F_Coun)+Par(F_Coun+1).*Par_Err(F_Coun+1))./C);
         Fi_Err = sqrt((Par_Err(F_Coun+1) + Par(F_Coun+1).*Par_Err(F_Coun)./Par(F_Coun))./((Par(F_Coun).^2+Par(F_Coun+1).^2)./Par(F_Coun))); 
         F_Line = [I,J,Har(I,1).*J,1./(Har(I,1).*J),Par(F_Coun),Par_Err(F_Coun),Par(F_Coun+1),Par_Err(F_Coun+1),C,C_Err,Fi,Fi_Err];
         fprintf(Fid,'%d  %d  %12.6f %12.6f   %10.5f %8.5f  %10.5f %8.5f  %10.5f %8.5f  %10.5f %8.5f\n',F_Line');
         F_Coun = F_Coun + 1;
      end
   end

   fprintf(Fid,'\n');
   fprintf(Fid,'%s \n','Linear Terms');
   fprintf(Fid,'%s \n','#D      Coef     Err');
   for I=0:1:Deg,
      F_Coun = F_Coun + 1;
      F_Line = [I, Par(F_Coun), Par_Err(F_Coun)];
      fprintf(Fid,'%d  %10.5f %8.5f \n',F_Line);
   end
   fprintf(Fid,'\n');
   fprintf(Fid,'\n');
   fprintf(Fid,'%s \n','  Fit quality:');
   fprintf(Fid,'%s ',    ['No. Deg. of Fredom : ']);
   fprintf(Fid,' %d \n',Freedom);
   fprintf(Fid,'%s',['              Chi2 : ']);
   fprintf(Fid,' %10.4f \n',Chi2);
   fprintf(Fid,'%s',['      Reduced Chi2 : ']);
   fprintf(Fid,' %10.4f \n',Chi2/Freedom);
   fprintf(Fid,'%s', ['   Cumulative Chi2 : ']);
   fprintf(Fid,' %6.4f \n',chi2cdf(Chi2,Freedom));

   fclose(Fid);
end



% calculating amplitude and phase
Nab = 2.*sum(Har(:,2));
for I=1:2:Nab-1,
   A   = Par(I);
   B   = Par(I+1);
   DA  = Par_Err(I);
   DB  = Par_Err(I+1);
   % calculate amplitude
   C   = sqrt(A.^2+B.^2);
   DC  = sqrt(((A.*DA).^2+(B.*DB).^2)./(A.^2+B.^2));  
   % calculate phase
   Ph  = atan2(B,A);
   DPh = sqrt((A.*DB).^2+(B.*DA).^2)./(A.^2+B.^2);
   % convert phase from radian to fraction
   Ph  = Ph./(2.*pi);
   DPh = DPh./(2.*pi);

   Par1(I,1)   = C;
   Par1(I,2)   = DC;
   Par1(I+1,1) = Ph;
   Par1(I+1,2) = DPh;
end
