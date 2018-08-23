function [StD,Theta,Bias,Z0,Ac]=jackknife(Data,Function,varargin)
% Given an estimator, calculate the Jacknife StD.
% Package: Util.stat
% Description: Given an estimator (given by a function), calculate the
%              Jackknife StD and the first order Quenouille-Tukey
%              Jacknife bias for this estimator.
%              Notes: - the estimator should be continues.
%                     - if Bias/StD~<0.25 then the bias is probably not
%                       an issue.
%                     - The bias estimate is not reliable if the estimator
%                       is an unsmooth statistic (e.g., median).
% Input  : - vector of data set.
%          - Estimator function name,
%            Fun(Data,Additional_Parameters)
%          - Additional parameters of 'Fun'.
% Output : - Jackknife StD.
%          - Vector of jackknife estimator.
%            That could be used for plotting the estimator distribution,
%            or to find percentile confidence interval using err_cl.m
%          - Jackknife Bias.
%          - bc_a bias correction (Z0).
%          - bc_a acceleration (Ac).
%            The bc_a bias coorection and acceleration can be used
%            to estimate the bias corrected and accelerated
%            confiedence interval.
%            For example: if you want the CI for the 0.9545 CL.
%                         use the function [A1,A2]=bc_a(Z0,Ac,0.9545),
%                         to estimate "corrected" percentile,
%                         and then CI=err_cl(Theta,0.9545);
% Reference : Efron B., 1982, in: The Jackknife, the bootstrap and other
%             resampling plans.
%             Efron B. & Tibshirani, R.J., 1993,
%             in: An Introduction to the Bootstrap
% see also: bootstrap_std.m, bc_a.m, err_cl.m
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Sep 2002
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [StD,Bias]=Util.stat.jackknife(Data,'mean',1);
%          % then UnBiased_Theta = Theta - Bias
% Reliable: 1
%------------------------------------------------------------------------------


N = size(Data,1);
Nsim = N;
Next = length(varargin);
STR = [];
for I=1:1:Next
   STR = [STR,sprintf(',%e',varargin{I})];
end


Theta = zeros(N,1);
for I=1:1:N
   J = find([1:1:N].'~=I);
   Theta(I) = eval([Function,'(Data(J,1)',STR,')']);

end

MeanTheta = mean(Theta);

% No jacknife
ThetaN    = eval([Function,'(Data',STR,')']);



StD  = sqrt( ((N-1)./N) .* sum((Theta - MeanTheta).^2) );
Bias = (N - 1).*(MeanTheta - ThetaN);


%--- calculate the bias correction and acceleration ---
% Efron & Tibshirani, p. 186
P  = length(find(Theta<ThetaN))./Nsim;
Z0 = norminv(P,0,1);   % bias correction
% acceleration
Ac = sum((MeanTheta - Theta).^3)./(6.*(sum((MeanTheta - Theta).^2)).^(1.5));
%--- Use bc_a.m to estimate bias corrected and accelerated confiedence interval
