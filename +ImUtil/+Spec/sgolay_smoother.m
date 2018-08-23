function [XI,YISG]=sgolay_smoother(X,Y)
%------------------------------------------------------------------------------
% sgolay_smoother function                                              FitFun
% Description: Smooth an unevenly spaced function [X,Y] in which some of
%              the Y values may be NaNs (e.g., gaps in the data), using
%                a combination of interpolation and Savitzky-Golay
%                (polynomial) filtering.
%                Furthermore, the function attempts to find the largest
%                window size possible such that the residuals do not have
%                a trend.
% Input  : - Column vector of X values (unevnly spaced).
%          - Column vector of Y values correspinds to the X values.
%            Some of the values maybe NaNs.
%------------------------------------------------------------------------------


InPar.InterpMethod = 'linear';
InPar.Order        = 2;
InPar.Window       = 301;
InPar.GapsInterp   = 'poly';
InPar.PreMedFilt   = []; 5;

% sort by X
if (issorted(X))
   % do nothing
else
   [X,SI] = sort(X);
   Y      = Y(SI);
end

% check if X is evenly spaced
if (range(diff(X))>0)
   % non evenly spaced
   % interpolate
   
   Step = min(diff(X));
   XI    = (min(X):Step:max(X)).';
   YI    = interpolate(X,Y,XI,InPar.InterpMethod);

else
   % evenly space function
   XI    = X;
   YI    = Y;
end


if (~isempty(InPar.PreMedFilt))
   % pref SG median filtering for CR removal etc.
   YI = medfilt1(YI,InPar.PreMedFilt);
end


% interpolate over gaps
switch lower(InPar.GapsInterp)
 case 'interp'
    Inan       = find(isnan(YI));
    Inn        = find(~isnan(YI));
    Ynan       = interp1(XI(Inn),YI(Inn),XI(Inan));  
    YISG(Inan) = Ynan;
 case 'poly'
    Inan       = find(isnan(YI));
    Inrange    = find(isnan(medfilt1(YI,501)) & ~isnan(YI));
    plot(XI(Inrange),YI(Inrange),'r--')

    P = polyfit(XI(Inrange),YI(Inrange),1);
    YI(Inan) = polyval(P,XI(Inan));
 otherwise
    error('Unknown GapsInterp option');
end



YISG = YI;
YISG = sgolayfilt(YI,InPar.Order,InPar.Window);

