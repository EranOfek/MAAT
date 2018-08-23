function InLC=sf_interp(LC,BinSize,Time,InterpMethod,CCF)
%--------------------------------------------------------------------------
% sf_interp function                                            timeseries
% Description: Interpolate a time series, and estimate the errors due
%              to the interpolation using the structure function
%              of the time series. The error in each interpolated
%              point is calculated by adding in quadrature the error
%              of the neastest point with the amplitude of the
%              stracture function at the the lag equal to the
%              difference between the interpolated point and the
%              nearest point.
% Input  : - Observations [Time, Mag, Err].
%          - Structure function binning interval.
%          - Vector of times for which to interpolate the LC.
%          - Interpolation method, default is 'linear'.
%          - Optional CCF matrix (see ccf.m) with the
%            structure function in columns [1 6 7].
%            If not given then calculate the structure function.
% Output : - Interpolated LC, [Time, Mag, Err].
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Dec 2002
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference : Ofek & Maoz 2003 ApJ 594, 101
% Reliable: 2
%--------------------------------------------------------------------------

ColT = 1;
ColM = 2;
ColE = 3;

CalcCCF = 'n';

if (nargin==3),
   InterpMethod = 'linear';
   CalcCCF      = 'y';
elseif (nargin==4),
   CalcCCF      = 'y';
elseif (nargin==5),
   CalcCCF      = 'n';
else
   error('Illegal number of input arguments');
end


switch CalcCCF
 case 'y'
    CCF  = ccf(LC,LC,BinSize,'normal');
    SF   = CCF(:,[1 6 7]);
 case 'n'
    % allready given
    SF   = CCF(:,[1 6 7]);
 otherwise
    error('Unknown CalcCCF option');
end

Ig = ~isnan(SF(:,2));
SF = SF(Ig,:);
Ig = SF(:,2)<0;
SF(Ig,2) = 0;

Mag = interp1(LC(:,ColT),LC(:,ColM),Time);

Nt   = length(Time);
InLC = zeros(Nt,3);
InLC(:,1) = Time;
InLC(:,2) = Mag;
for I=1:1:Nt,
   DeltaT     = Time(I) - LC(:,ColT);
   [Min,MinI] = min(abs(DeltaT));
   %--- interpolate SF ---

   Amp2  = interp1(SF(:,1),SF(:,2),Min);
   Err   = sqrt(LC(MinI,ColE).^2 + Amp2);

   InLC(I,3) = Err;
end



