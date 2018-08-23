function [Peaks,Y,Thresh]=matched_filter(X,Filter,Subtract,Thresh)
% Matched filter for 1-D time series
% Package: timeseries
% Description: Run a matched filter (North filter) for a 1-D series and
%              return the possible peaks (signal) above a given threshold.
% Input  : - 1-D evenly spaced series.
%          - Filter. This is either a vector of filter, or a cell
%            array in which the first element is a function hendle
%            to generate the filter, and the other elements are input
%            arguments to pass to the filter function
%            (e.g., {'fun_gauss',[1 5 2],(1:1:10)',0}).
%          - Subtract a background prior to filtering.
%            This is a cell array of extra parameters to pass to
%            the subtract_back1d.m program (e.g., {'mean'}).
%            Default is {'none'}. If empty, then use default.
%          - Threshold above to select peaks. Default is 3.
%            If negative number is provided (e.g., -5), then will
%            estimate the sigma using the 68th percentile, and use
%            the absolute value of the threshold as the number of sigmas.
%            (e.g., -5 is 5 sigma).
% Output : - Matrix of values above threshold:
%            [X index, Y at index, X at index, X after background removal].
%          - The matched filter series (Y).
%          - Actual threshold used.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=ones(10000,1); X(5001:5005)=[1.5 2 2.2 2 1.5]';
%          X=X+randn(10000,1).*0.1; F=[1.5 2 2.2 2 1.5]';
%          [Peaks,Y]=timeseries.matched_filter(X,F);
%          [Peaks,Y,Thresh]=timeseries.matched_filter(X,F,{'median'},-5);
% Reliable: 2
%--------------------------------------------------------------------------

Def.Subtract = {};
Def.Thresh   = 3;
if (nargin==2),
    Subtract  = Def.Subtract;
    Thresh    = Def.Thresh;
elseif (nargin==3),
    Thresh    = Def.Thresh;
elseif (nargin==4),
    % do nothing
else
    error('Illegal number of input arguments');
end

if (iscell(Filter)),
    F = feval(Filter{1},Filter{2:end});
else
    F = Filter;
end
% normalize filter
F = F./sum(F);
Fhs = round((length(F)-1).*0.5);

if (isempty(Subtract)),
    Subtract = {'none'};
end
Xbs = timeseries.subtract_back1d(X,Subtract{:});


Y  = abs(ifft(fft(Xbs).*fft(F,length(Xbs))));

if (Thresh<0),
    Thresh = quantile(Y,0.6827).*abs(Thresh);
end

Ip = find(Y>Thresh);
Ip     = Ip((Ip-Fhs)>1);  % remove peaks near edges
Peaks=[Ip, Y(Ip), X(Ip-Fhs), Xbs(Ip-Fhs)];
