function D=wmedian(varargin)
% Weighted median for a vector.
% Package: Util.stat
% Description: Weighted median for a vector.
%              Calculates the weighted median of a vector
%              given the error on each value in the vector.
% Input  : * Either two arguments Data,Error or a two column matrix
%            of [Data,Error].
% Output : - Weighted median
% Reference: https://en.wikipedia.org/wiki/Weighted_median
% License: GNU general public license version 3
% Tested : Matlab R2015a
%     By : Eran O. Ofek                    Jun 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: D=wmedian(V,E);
% Reliable: 2
%--------------------------------------------------------------------------

if (numel(varargin)==1)
    Val = varargin{1}(:,1);
    Err = varargin{1}(:,2);
else
    Val = varargin{1};
    Err = varargin{2};
end

[SV,SI] = sort(Val);
CSW     = cumsum(1./Err(SI).^2);  % cumsum of the weights
if (any(isinf(CSW)))
    D = NaN;
else
    CSWnorm = CSW./CSW(end);
    D = interp1(CSWnorm,SV,0.5,'linear');
end

    