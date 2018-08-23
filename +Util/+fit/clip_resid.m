function [I,NR,Flag]=clip_resid(Resid,varargin)
% Clip residuals using various methods.
% Package: Util.fit
% Description: Clip residuals using various methods including sigma
%              clipping min/max, etc.
% Input  : - Vector of residuals to clip.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            The following keywords are available:
%            'Method'  - The method of clipping. Available methods are:
%                        'StD'      - Sigma clipping, where sigma is
%                                     estimated based on the StD (default).
%                        'StdP'     - Sigma clipping, where sigma is
%                                     estimated based on the 68 percentile.
%                                     In this case the std is not necesserly
%                                     symmetric.
%                        'Constant' - Lower/upper constant above/below the
%                                     mean.
%                        'MinMax'   - Reject N lowest and highest values.
%                        'Perc'     - Reject lower and upper percentiles.
%                                     This method ignores the 'Mean' option,
%                                     and the rejection is based on percentile.
%            'Mean'    - Method by which to calculate the "mean" of the
%                        sample. Available methods are:
%                        'mean'   - mean of the sample.
%                        'median' - median of the sample (default).
%            'Clip'    - Two elements vector containing the lower and upper
%                        values for the sigma clipping.
%                        This is [Lower, Upper] number of sigmas (positive)
%                        below/above the mean.
%                        For Method = 'StD'  or 'StdP' this is the number
%                        of sigma below/above the mean to reject.
%                        For 'Constant' this is the number below/above
%                        the mean to reject.
%                        For 'MinMax' this is the number of lowest/highest
%                        residuals to reject.
%                        For 'Perc' this is the lower/upper percentile
%                        fraction to reject. For example [0.02 0.02],
%                        will reject the upper and lower 2%.
%                        Default is always set to no rejection
%                        (i.e., [Inf Inf] or [0 0] depending on Method.
%                        Empty matrix will use default.
%            'StdZ'    - Add epsilon to std {'y'|'n'}, default is 'y'.
%                        This is useful in case std is zero.
% Output : - Index of good residuals.
%          - Vector of good residuals.
%          - For each residual, a flag indicating if residual is good (1) or
%            outside cliped region(0).
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    May 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: I=clip_resid(randn(100,1)+5,'Method','stdp','clip',[3 3]);
% Reliable: 2
%--------------------------------------------------------------------------
ONE_SIGMA = 0.6827;

DefV.Method   = 'StD';
DefV.Mean     = 'median';
DefV.Clip     = [];
DefV.StdZ     = 'y';
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});   % ignore additional keywords...
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

switch lower(InPar.StdZ)
 case 'y'
    AddStd0 = eps;
 case 'n'
    AddStd0 = 0;
 otherwise
    error('Unknown Std0 option');
end

switch lower(InPar.Mean)
 case 'mean'
    Mean = nanmean(Resid);
 case 'median'
    Mean = nanmedian(Resid);
 otherwise
    error('Unknown Mean option');
end

if (isempty(InPar.Clip)),
   switch lower(InPar.Method)
    case {'std','stdp','constant'}
       InPar.Clip = [Inf Inf];
    case {'minmax','perc'}
       InPar.Clip = [0 0];
    otherwise
       error('Unknown Method option');
   end
end



switch lower(InPar.Method)
 case 'std'
    Std   = nanstd(Resid) + AddStd0;
    Lower = Mean - InPar.Clip(1).*Std;
    Upper = Mean + InPar.Clip(2).*Std;
 case 'stdp'
    Std   = err_cl(Resid,ONE_SIGMA) - Mean;  
    Lower = Mean + InPar.Clip(1).*Std(1) - AddStd0;
    Upper = Mean + InPar.Clip(2).*Std(2) + AddStd0;
 case 'constant'
    Lower = Mean - InPar.Clip(1);
    Upper = Mean + InPar.Clip(2);
 case 'minmax'
    if (InPar.Clip(1)==0),
       Lower = Inf;
    end
    if (InPar.Clip(2)==0),
       Upper = Inf;
    end
    if (InPar.Clip(1)==0 && InPar.Clip(2)==0),
       % do nothing
    else
       [SortedResid] = sort(Resid);
       Lower = SortedResid(InPar.Clip(1))+eps;
       Upper = SortedResid(end-InPar.Clip(2)+1)-eps;
    end
 case 'perc'
    Lower = prctile(Resid,InPar.Clip(1).*100);
    Upper = prctile(Resid,(1-InPar.Clip(2)).*100);
 otherwise
    error('Unknwon Method option');
end


Flag = Resid>Lower & Resid<Upper;
I = find(Flag);
if (nargout>1),
   NR = Resid(I);
end

 


