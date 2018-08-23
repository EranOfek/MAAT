function [StD,Theta,Bias,Z0,Ac]=bootstrap_std(Data,Function,Nsim,varargin)
%--------------------------------------------------------------------------
% bootstrap_std function                                         AstroStat
% Description: Given an estimator (given by a function) calculate the
%              Bootstrap StD for this estimator.
% Input  : - Matrix of data set, row per measurment.
%          - Estimator function name,
%            Fun(Data,Additional_Parameters)
%          - Number of Bootstrap simulations.
%          - Additional parameters of 'Fun'.
% Output : - Bootstrap StD.
%          - Vector of bootstrap estimator.
%            That could be used for plotting the estimator distribution,
%            or to find percentile confidence interval using err_cl.m
%          - Bootstrap estimate of the bias.
%            Subtract the bias from the estimator
%            to get the corrected estimator.
%          - bc_a bias correction (Z0).
%          - bc_a acceleration (Ac).
%            The bc_a bias coorection and acceleration can be used
%            to estimate the bias corrected and accelerated
%            confiedence interval.
%            For example: if you want the CI for the 0.9545 CL.
%                         use the function [A1,A2]=bc_a(Z0,Ac,0.9545),
%                         to estimate "corrected" percentile,
%                         and then CI=err_cl(Theta,0.9545);
% Reference : Efron B. & Tibshirani, R.J., 1993,
%             in: An Introduction to the Bootstrap
% Example: [StD,Theta]=bootstrap_std(Data,'mean',1000,1)
%          then UnBiased_Theta = Theta - Bias
% see also: jackknife.m, bc_a.m, err_cl.m
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 2002
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [StD,Theta,Bias,Z0,Ac]=bootstrap_std(randn(100,1),@mean,1000);
% Reliable: 1
%--------------------------------------------------------------------------


N = size(Data,1);
Next = length(varargin);
STR = [];
for I=1:1:Next,
   STR = [STR,sprintf(',%e',varargin{I})];
end

ParIndex = 1;

Theta = zeros(Nsim,1);
for I=1:1:Nsim,
   Rand    = rand(N,1);
   NewInd  = ceil(Rand.*N);

   NewData = Data(NewInd,:);

   %   Theta(I) = eval([Function,'(NewData(:,:)',STR,')']);
   PP = feval(Function,NewData,varargin{:});
   Theta(I) = PP(ParIndex); %feval(Function,NewData,varargin{:});

end

%ThetaN = eval([Function,'(Data(:,:)',STR,')']);
%PP = feval(Function,NewData(:,:),varargin{:});
PP = feval(Function,Data(:,:),varargin{:});
ThetaN = PP(ParIndex); %feval(Function,NewData(:,:),varargin{:});

MeanTheta = mean(Theta);

StD       = sqrt(sum((Theta - MeanTheta).^2)./(Nsim-1));

Bias      = MeanTheta - ThetaN;


%--- calculate the bias correction and acceleration ---
% Efron & Tibshirani, p. 186
P  = length(find(Theta<ThetaN))./Nsim;
Z0 = norminv(P,0,1);   % bias correction
% acceleration
Ac = sum((MeanTheta - Theta).^3)./(6.*(sum((MeanTheta - Theta).^2)).^(1.5));
%--- Use bc_a.m to estimate bias corrected and accelerated confiedence interval
