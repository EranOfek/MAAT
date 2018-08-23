function [NewData,Nbounds]=error2ensemble(Data,Errors,Fun,Bounds)
%--------------------------------------------------------------------------
% error2ensemble function                                        AstroStat
% Description: Generate a realization of data points given the data
%              probability distribution function.
% Input  : - Column vector of the expectency value for each data point.
%          - One or two column matrix containing the lower (first column)
%            and upper (second column) error for each data point
%            realization. If single column is given then the lower and
%            upper errors are equal.
%            If empty matrix (i.e., []), then the distribution function
%            is given numerically.
%          - Function to use in order to generate the realizations.
%            This can be a mumerical function [X, P], where X=0
%            corresponds to the expectency value of the distribution.
%            Alternatively this can be a string containing a pre-defined
%            function: {'normal'|'lognormal'|'poisson'|'flat'},
%            default is 'normal'.
%            'normal'  use the lower and upper errors as the lower and
%                      upper 1\sigma errors.
%            'flat'    use the errors as the bounderies of the flat
%                      distribution.
%            'poisson' use the first column in errors as lambda
%                      (the second column is ignored).
%          - Two columns matrix with the [minimum, maximum] bounds
%            of the minumum and maximum allowed in each realization.
%            Default is [-Inf Inf].
%            If a data point with value below or above these bounds
%            is found, then it is being replaced by the bound value.
% Output : - Column vector of data realizations.
%          - Number of boundary replacments.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    May 1998
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [NewData,Nbounds]=error2ensemble(ones(1000,1),[1 5]);
% Reliable: 2
%--------------------------------------------------------------------------
DefFun       = 'normal';
DefBounds    = [-Inf Inf];

if (nargin==2),
   Fun      = DefFun;
   Bounds   = DefBounds;
elseif (nargin==3),
   Bounds   = DefBounds;
elseif (nargin==4),
   % no default
else
   error('Illegal number of input arguments');
end

if (size(Errors,2)==1),
   Errors = [Errors, Errors];
end

N = length(Data);

if (size(Errors,1)==1),
   Errors = ones(N,1)*Errors;
end

if (size(Bounds,1)==1),
   Bounds = ones(N,1)*Bounds;
end

if (ischar(Fun)==1),
   switch lower(Fun)
    case 'normal'
       Rand = randn(N,1);
       I1   = find(Rand<0);
       I2   = find(Rand>=0);
       NewData = Data;
       if (isempty(I1)==0),
          NewData(I1) = NewData(I1) + Rand(I1).*Errors(I1,1);
       end
       if (isempty(I2)==0),
          NewData(I2) = NewData(I2) + Rand(I2).*Errors(I2,2);
       end
    case 'lognormal'
       error('lognormal not implemented yet');

    case 'poisson'
       NewData   = poissrnd(Errors(:,1),N,1);

    case 'flat'
       NewData = rand(N,1).*(Errors(:,2)-Errors(:,1)) + Errors(:,1);

    otherwise
       error('Unknown Fun option');
   end
else
   % numerical function
   error('numerical function not implemented yet');
end


I1 = find(NewData<Bounds(:,1));
I2 = find(NewData>Bounds(:,2));
NewData(I1) = Bounds(I1,1);
NewData(I2) = Bounds(I2,2);
Nbounds     = length(I1) + length(I2);

