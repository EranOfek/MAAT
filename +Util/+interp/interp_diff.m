function [Yi,SPoly,H,X0]=interp_diff(X,Y,Xi,Deg,Check)
% Interpolation based on 4th order Stirling formula
% Package: Util.interp
% Description: Interpolation of equally spaced data using 
%              high-order differences.
% Input  : - Equally spaced and (asendingly) sorted X.
%          - Y
%          - X values for which to interpolate.
%          - Degree of differences, default is 4.
%          - Check if X is equally spaced {true|false}, default is false.
% Output : - Interpolated Y values.
%            Return NaN in case that extrapolation is needed.
%          - Stirling interpolation polynomial (only for 4th deg).
%            Each line contains a4, a3, a2, a1, a0 coeff.
%            then: f_p=(((a4*p+a3)*p+a2)*p+a1)*p+a0 ,
%            where p is the interpolation factor: p=(X-X0)./H 
%            The polynomial is suitable for use in the range -0.5<p<0.5,
%            and it may be adequate in the range -2<p<2. 
%          - The tabulation interval H.
%          - X0 of the interpolation for each point/polynomial.
% See also: find_local_extramum.m, find_local_zeros.m, interp_diff_ang.m
% Reference: Seidelmann 1992, Explanatory Supp. to the Astron. Almanac
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=[-5:1:5].';  Y=(X-0.2).^2;
%          [Yi,SPoly,H,X0]=Util.interp.interp_diff(X,Y,0.2);
% Reliable: 1
%--------------------------------------------------------------------------
DefDeg     = 4;
DefCheck   = false;
if (nargin==3)
   Deg    = DefDeg;
   Check  = DefCheck;
elseif (nargin==4)
   Check  = DefCheck;
elseif (nargin==5)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (size(X,1)>=size(X,2))
   % do nothing
else
   X = X.';
   Y = Y.';
end

Deg   = 4;
Check = false;

N  = length(Y);
if (isempty(X))
    X = (1:1:N).';
end

if (Check)
    if (length(find(diff(X,2,1)~=0))>0 || (X(2)<X(1)))
       error('X must be equally spaced');
    end
end

Ni = length(Xi);
Yi = zeros(Ni,1);
X0 = zeros(Ni,1);
H  = X(2) - X(1);

%FlagOutOfRange     = (Xi-min(X))./H<2 | (Xi-max(X))./H>(N-3);
%Yi(FlagOutOfRange) = NaN;
%X0(FlagOutOfRange) = NaN;

Diff1 = [diff(Y,1); NaN];
Diff2 = [NaN; diff(Y,2); NaN];
Diff3 = [NaN; diff(Y,3); NaN; NaN];
Diff4 = [NaN; NaN; diff(Y,4); NaN; NaN];

Ix    = zeros(Ni,1) + 2;
for Ii=1:1:Ni
    if (~isnan(Yi(Ii)))
        Ix(Ii) = find([X;Inf]>=Xi(Ii),1)-1;
        %X0(Ii) = X(Ix);
        %Y0(Ii) = Y(Ix);
    end
end
Ix(Ix<2) = 2;
Ix(Ix>(N-1)) = N-1;
X0 = X(Ix);

p = (Xi - X0)./H;
B2 = p.*(p-1).*0.25;
B3 = B2.*(p-0.5).*2./3; %% p.*(p-1).*(p-0.5)./6;
B4 = (p+1).*p.*(p-1).*(p-2)./48;

Yi = Y(Ix) + p.*Diff1(Ix) + B2.*(Diff2(Ix)+Diff2(Ix+1)) + B3.*Diff3(Ix) + B4.*(Diff4(Ix)+Diff4(Ix+1));

if (nargout>1)
    SPoly = zeros(Ni,Deg+1);    % [A4 A3 A2 A1 A0]
    SPoly(:,5) = Y(Ix);                                      % A0
    SPoly(:,2) = (Diff3(Ix) + Diff3(Ix-1))./12;              % A3
    SPoly(:,1) = Diff4(Ix)./24;                              % A4
    SPoly(:,3) = Diff2(Ix).*0.5 - SPoly(:,1);                % A2
    SPoly(:,4) = (Diff1(Ix)+Diff1(Ix-1)).*0.5 - SPoly(:,2);  % A1
end




