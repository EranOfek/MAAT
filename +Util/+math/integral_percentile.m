function [LeftLimit,RightLimit,TotInt]=integral_percentile(X,Y,Per,Type,Method,Norm)
% Given a tabulate function find limits that contains percentile
% Package: Util.math
% Description: Given a numerically tabulated function (X,Y) and a
%              percentile P, find the limits (a,b) such that
%              int_{a}^{b}(Y(X)dx)=P (i.e., the integral of Y(X) from a
%              to b equal P.
% Input  : - X
%          - Y
%          - Vector of percentile, for each to calculate the limits
%            of the integral which contains that percentiles of the
%            total integral.
%          - Percentile type:
%            'c'  - central (default). The limits are such that the
%                   integarl from min(X) to the left limit is equal
%                   to the integral from the right limit to max(X).
%            'l'  - The limits are such that the left limit is 0.
%            'r'  - The limits are such that the right limit is 0.
%          - Method:
%            'Interp'  - interpolate the integral of the tabulated
%                        function.
%                        This option works well for smooth functions
%                        (default).
%            'First'   - Do not interpolate the integral of the
%                        tabulated function and select the first
%                        point that satisfy the integral.
%                        This is recomended only if the function
%                        tabulation is dense.
%          - Normalize the tabulated function such its integral equal
%            to 1 {'y' | 'n'}. Default is 'y'.
% Output : - Vector of left limits.
%          - Vector of right limits.
%          - Total integral.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Dec 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%-------------------------------------------------------------------------
InterpMethod = 'linear';

Def.Type   = 'c';
Def.Method = 'Interp';
Def.Norm   = 'y';

if (nargin==3)
   Type   = Def.Type;
   Method = Def.Method;
   Norm   = Def.Norm;
elseif (nargin==4)
   Method = Def.Method;
   Norm   = Def.Norm;
elseif (nargin==5)
   Norm   = Def.Norm;
elseif (nargin==6)
   % do nothing
else
   error('Illegal number of input arguments');
end


CumInt    = cumtrapz(X,Y);
TotInt    = max(CumInt);
PerCumInt = CumInt./TotInt;

switch lower(Norm)
 case 'y'
    % do nothing
 case 'n'
    PerCumInt = CumInt;
 otherwise
    error('Unknown Norm option');
end

CumInt

switch lower(Type)
 case 'c'
    switch lower(Method)
     case 'interp'
        LeftLimit   = interp1(PerCumInt,X,    0.5.*(1-Per), InterpMethod);
        RightLimit  = interp1(PerCumInt,X,1 - 0.5.*(1-Per), InterpMethod);
     case 'first'
        for Ip=1:1:length(Per)
           I = find( (PerCumInt - 0.5.*(1-Per(Ip))) > 0);
           LeftLimit(Ip,1) = X(I(end));
           I = find( (PerCumInt - (1-0.5.*(1-Per(Ip)))) < 0);
           RightLimit(Ip,1) = X(I(1));
        end
     otherwise
        error('Unknown Method Option');
    end
 case 'l'
    RightLimit = interp1(PerCumInt,X,Per, InterpMethod);
    LeftLimit  = zeros(size(RightLimit));
 case 'r'
    LeftLimit  = interp1(PerCumInt,X,1-Per, InterpMethod);
    RightLimit = ones(size(LeftLimit)).*max(X);
 otherwise
    error('Unknown Type Option');
end
