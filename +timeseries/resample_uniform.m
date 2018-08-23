function [IntVal,NewT]=resample_uniform(T,Val,varargin)
% Uniform resampling of a non-evenly spaced time series.
% Package: timeseries
% Description: Uniform resampling of a non-evenly spaced time series.
% Input  : - A column vector of times.
%          - A matrix of values to resample, one column per property.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'DeltaT' - If a scalar than this is the time step (in the tim
%                       units) of the equally spaced resampled time series.
%                       If a vector than this is the vector of time at
%                       which to resample the light curve (not necesserly
%                       evenly spaced). Default is 1.
%            'InterpMethod' - A cell array of interpolation methods with
%                       which to resample the light curve. A method per
%                       column in the second input argument. If a single
%                       method is provided than it will be applied to all
%                       columns.
%                       Default is 'linear'.
% Output : - A matrix of resampled light curves (one per column).
%          - A vector of the resampled times.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [IntVal,NewT]=timeseries.resample_uniform(rand(100,1).*100,rand(100,1))
% Reliable: 2
%--------------------------------------------------------------------------


DefV.DeltaT               = 1;
DefV.InterpMethod         = 'linear';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (~iscell(InPar.InterpMethod))
    InPar.InterpMethod = {InPar.InterpMethod};
end
NInt = numel(InPar.InterpMethod);

Ncol = size(Val,2);

% sort
[~,SI] = sort(T);
Val    = Val(SI,:);
T      = T(SI);

% set new resampling grid
if (numel(InPar.DeltaT)>1)
    NewT = InPar.DeltaT(:);
else
    NewT = (T(1):InPar.DeltaT:T(end)).';
end
    


IntVal = zeros(numel(NewT),Ncol);
for Icol=1:1:Ncol
    
    % resample
    Ic = min(NInt,Icol);
    switch lower(InPar.InterpMethod{Ic})
        case {'linear','nearest','next','previous','spline','pchip','cubic','v5cubic'}
            % use interp1
            IntVal(:,Icol) = interp1(T,Val(:,Icol),NewT,InPar.InterpMethod{Ic});
        otherwise
            error('Uknown InterpMethod option');
    end
end

