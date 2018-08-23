function [MaxX,MaxY]=find_peak_center(X,Y,varargin)
%--------------------------------------------------------------------------
% find_peak_center function                                        General
% Description:
% Input  : - Evenly spaced X position.
%            If empty then will use 1 to length of Y.
%          - Y values.
%          * Arbitrary number of ...,key,val,... input arguments.
%            The following keywords are available:
%            'Method' - One of the following fitting method:
%                       'max'    - maximum
%                       'wmean'  - Mean weighted by inverse of counts (Y).
%                       'wmeanbs'- Mean weighted by inverse of counts (Y)
%                                  after median background subtraction.
%                       'str'    - Striling interpolation.
%                                  Return all peaks.
%                       'strh'   - Striling interpolation.
%                                  Return highest peak only.
%                       'strn'   - Striling interpolation.
%                                  Return peak nearest to center (default).
%            'Ind' - Indices of X/Y vectors to use in search.
%                    If empty use all. Default is empty.
% Output : - Max X value. Return NaN if no peak found.
%          - Max Y value. Return NaN if no peak found.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=(1:1:10).'; Y=ones(size(X)); Y(4:6)=[2 3 2.1];
%          [MaxX,MaxY]=find_peak_center(X,Y,'Method','strh');
% Reliable: 2
%--------------------------------------------------------------------------



DefV.Method      = 'strn';
DefV.Ind         = [];
%DefV.FitFun      = fun_gauss([1 5 1],(1:1:10).');
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

N = numel(Y);
if (isempty(X)),
    X = (1:1:N).';
end

if (isempty(InPar.Ind)),
    FI = true(N,1);
else
    FI = InPar.Ind;
end

X = X(FI);
Y = Y(FI);
switch lower(InPar.Method)
    case 'max'
        [MaxY,MaxInd] = max(Y);
        MaxX = X(MaxInd);
    case 'wmean'
        MaxX = sum(X.*Y)./sum(Y);
        MaxY = max(Y);
    case 'wmeanbs'
        Ybs  = Y - median(Y);
        MaxX = sum(X.*Ybs)./sum(Ybs);
        MaxY = max(Y);
    case {'str','strh','strn'}
        Extram=find_local_extramum(X,Y);
        %[X, Y, 2nd derivative d^2Y/dX^2]
        % look for local maxima
        Extram = Extram(Extram(:,3)<0,:);
        MaxX   = Extram(:,1);
        MaxY   = Extram(:,2);
        switch lower(InPar.Method)
            case 'strh'
                [~,Ihp] = max(MaxY);
                MaxX = MaxX(Ihp);
                MaxY = MaxY(Ihp);
            case 'strn'
                % select peak nearest to center
                [~,Inp] = min(abs(MaxX - mean(X)));
                MaxX = MaxX(Inp);
                MaxY = MaxY(Inp);
            otherwise
                % do nothing
        end
     
    otherwise
        error('Unknown Method option');
end
if (isempty(MaxX)),
    MaxX = NaN;
    MaxY = NaN;
end