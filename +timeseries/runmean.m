function RunMean=runmean(Data,WFun,WFunPars,WS,TimeVec)
% Running mean of un-evenly spaced time series
% Package: timeseries
% Description: Calculate the runing mean of an unevenly spaced time
%              series with different weight functions and weight schemes.
% Input  : - Data matrix, [Time, Value, Error].
%            Error is optional, if not given then assumes
%            equal weights.
%          - Weight function type:
%            'f' : flat function:
%                  parameter is : i) total width
%            'g' : gaussian function (truncated at 3\sigma):
%                  parameter is : i) HWHM width
%            'c' : cosine-bell function:
%                  parameters are: i) total width
%                                 ii) fraction of flat part
%          - Weight function vector of parameters.
%          - Weighting scheme:
%            'f2' - use weight function squared only, default.
%            'f'  - use weight function only.
%            'wm' - use points-errors (weighted mean) only.
%            'wf' - use sum of squares of weight function and
%                   normalize points-errors so their mean is RelWeight.
%          - Time vector for which to calculate the running mean,
%            default is Data(:,1).
% Output : - Runing mean matrix [Time, Value].
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jan 2002
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: timeseries.runmean(rand(100,3),'c',[0.1 0.05],'wf');
% Reliable: 2
%------------------------------------------------------------------------------

RelWeight = 0.1;    % relative weight between errors and weight function
TCol = 1;
YCol = 2;
ECol = 3;
if (nargin==4),
   TimeVec = Data(:,TCol);
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end

DataSpan = max(Data(:,1)) - min(Data(:,1));
SizeData = size(Data);
if (SizeData(2)==2),
   Data = [Data, ones(SizeData(1),1)];
end


% build weight function
switch WFun
 case 'f' 
    % flat function
    X     = [-0.5.*WFunPars(1):0.01.*WFunPars(1):0.5.*WFunPars(1)].';
    Y     = ones(size(X));
   
 case 'g'
    % truncated gaussian
    SC    = 3; % sigma cutoff  
    HWHM  = WFunPars(1);
    Sigma = sqrt(HWHM.^2./(2.*log(2)));
    X     = [-SC.*Sigma:0.02.*Sigma:SC.*Sigma].';
    Y     = exp(-X.^2./(2.*Sigma.^2));
    HalfR = SC.*Sigma;  % half range
    
 case 'c' 
    % cosine bell
    Range        = [-0.5.*WFunPars(1):0.01.*WFunPars(1):0.5.*WFunPars(1)].';
    [Y,X] = timeseries.cosbell(WFunPars(2),Range);
    
 otherwise
    error('Unknown WFun type');
end

NormErr = RelWeight.*Data(:,ECol)./mean(Data(:,ECol));
InterpMethod = 'linear';
N  = length(TimeVec);
RM = zeros(N,1);
for I=1:1:N,
   if (isnan(TimeVec(I))==0),
      % points weight
      Wp = interp1(X+TimeVec(I),Y,Data(:,TCol),InterpMethod);
      
      switch WS
       case 'f2'
          % 'f2' - use weight function squared only, default.
          RM(I) = nansum(Wp.^2 .* Data(:,YCol))./nansum(Wp.^2);
          
       case 'f'
          % 'f'  - use weight function only.
          RM(I) = nansum(Wp .* Data(:,YCol))./nansum(Wp);
          
       case 'wm'
          % 'wm' - use points-errors (weighted mean) only.
          RM(I) = nansum(ceil(Wp) .* Data(:,YCol) .* (1./Data(:,ECol)).^2)./nansum(ceil(Wp) .* (1./Data(:,ECol)).^2);
          
       case 'wf'
          % 'wf' - use sum of squares of weight function and
          %        normalize points-errors so their mean is 1.
          RM(I) = nansum(Data(:,YCol) .* ((1./NormErr).^2 + Wp.^2))./nansum((1./NormErr).^2 + Wp.^2);
          
       otherwise
          error('Unknown Weights scheme type');
      end
   else
      RM(I) = NaN;
   end   
end

RunMean = [TimeVec, RM];
