function [XCS,InfoX]=xcorr_stretch(S1,S2,StretchVec,Back,varargin)
% Find the stretch and shift between two seriese using fft.
% Package: Util.xcorr
% Description: Find the stretch and shift between two seriese using
%              a stretch and fft approach.
%              In order to match the first vector to the second (reference)
%              vector one need to stretch the first vector and then shift
%              it by the amount found by this script.
% Input  : - Vector of first series (equaly spaced).
%            If this is a two column matrix, then the first column
%            is (equally spaced) X-position, and the second column
%            is the Y value.
%          - Second vector (reference).
%            See first input argument for details.
%          - Vector of stretches to test.
%            Default is [0.1:0.01:10]';
%          - Background subtraction method:
%            'none' - do not subtract background.
%            'mean' - subtract the mean.
%            'median' - subtract the median (default).
%            'medfilt' - subtract a median filter, which length is
%                        specified as the next input argument.
%          * Arbitrary number of input argument to pass to the
%            background subtraction algorithm.
% Output : - Matrix with the following columns:
%            [Stretch, Correlation, Shift].
%            Note that the stretch is operating prior to the shift.
%          - Structure containing information about te best stretch
%            and shift between the vectors. Fields are:
%            .MaxCorr   - Maximum correlation.
%            .X2toX1    - [Stretch, Shift] needed to transform
%                         the [X2 S2] series to the frame of
%                         the [X1 S1] series.
%            .X1toX2    - [Stretch, Shift] needed to transform
%                         the [X1 S1] series to the frame of
%                         the [X2 S2] series.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    August 2010
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Algorithm: Given a refernce series [X2, S2]
%            and a shifted/stretched series [X'2, S2] (AKA [X1, S1]).
%            Let the hidden relation (which we want to find)
%            between the two serieses be: X'2 = alpha*X2 + beta.
%            Let Xn2=[1:1:N2] and Xn1=[1:1:N1], and
%            Xn2 = alpha_n*X2 + beta_n ; Xn1 = alpha_1*X'2 + beta_1
%            [alpha_1, beta_1] = polyfit(X'2,Xn1);
%            [alpha_n, beta_n] = polyfit(X2,Xn2);
%            Next we have to find using the stretch and cross-correlate
%            method gamma and delta such:.
%            Xn2 = gamma*Xn1 + delta.
%            The solution is:
%            alpha = alpha_n/(gamma*alpha_1);
%            beta  = (beta_n + gamma*beta_1 + delta)/(gamma.*alpha_1);
% Example: X1 = [1:1:1000].'; X2 = X1.*1.5+100;
%          R1 = rand(1e3,1); R2 = rand(1e3,1);
%          R1(11:20) = 10; R1(111:120) = 15;
%          R2(201:220) = 9; R2(401:419) = 16;
%          [XCS,Info] = Util.xcorr.xcorr_stretch(R1,R2);
%          % or with X coordinates:
%          [XCS,Info] = Util.xcorr.xcorr_stretch([X1 R1],[X2 R2]);
%---------------------------------------------------------------------------
Def.StretchVec = [0.1:0.01:10].';
Def.Back = 'median';
if (nargin==2),
   StretchVec = Def.StretchVec;
   Back       = Def.Back;
elseif (nargin==3),
   Back       = Def.Back;
end

if (size(S1,2)>1),
   X1 = S1(:,1);
   S1 = S1(:,2);
   N1 = length(S1);
else
   N1 = length(S1);
   X1 = [1:1:N1]';
end

if (size(S2,2)>1),
   X2 = S2(:,1);
   S2 = S2(:,2);
   N2 = length(S2);
else
   N2 = length(S2);
   X2 = [1:1:N2]';
end

Xn1 = [1:1:N1].';
Xn2 = [1:1:N2].';

switch lower(Back)
 case 'none'
    % do nothing
 case 'median'
    S1 = S1 - median(S1);
    S2 = S2 - median(S2);
 case 'mean'
    S1 = S1 - mean(S1);
    S2 = S2 - mean(S2);
 case 'medfilt'
    S1 = S1 - medfilt(S1,varargin{:});
    S2 = S2 - medfilt(S2,varargin{:});
 otherwise
    error('Unknown Back option');
end

Nsv        = length(StretchVec);
XCS        = zeros(Nsv,3);
for Isv=1:1:Nsv,
   SS1 = interp1([1:1:N1].',S1,1./StretchVec(Isv).*[1:1:N1].');
   SS1(find(isnan(SS1))) = 0; % replace NaN with zeros

   [Lag,XC,Info] = xcorr_fft(SS1,S2,'none');

   XCS(Isv,:) = [StretchVec(Isv), Info.Corr, -Info.BestShift];
end
   
% find best shift and stretch
[Max,MaxInd]   = max(XCS(:,2));
InfoX.MaxCorr  = Max;
% stretch and shift in [1:1:N] system:
Gamma          = XCS(MaxInd,1);
Delta          = XCS(MaxInd,3);


% find transformation between original coordinate systems:
% Algorithm: Given a refernce series [X2, S2]
%            and a shifted/stretched series [X'2, S2] (AKA [X1, S1]).
%            Let the hidden relation (which we want to find)
%            between the two serieses be: X'2 = alpha*X2 + beta.
%            Let Xn2=[1:1:N2] and Xn1=[1:1:N1], and
%            Xn2 = alpha_n*X2 + beta_n ; Xn1 = alpha_1*X'2 + beta_1
%            [alpha_1, beta_1] = polyfit(X'2,Xn1);
%            [alpha_n, beta_n] = polyfit(X2,Xn2);
%            Next we have to find using the stretch and cross-correlate
%            method gamma and delta such:.
%            Xn2 = gamma*Xn1 + delta.
%            The solution is:
%            alpha = alpha_n/(gamma*alpha_1);
%            beta  = (beta_n + gamma*beta_1 + delta)/(gamma.*alpha_1);

Par1    = polyfit(X1,Xn1,1);
Alpha_1 = Par1(1);
Beta_1  = Par1(2);
Par2    = polyfit(X2,Xn2,1);
Alpha_n = Par2(1);
Beta_n  = Par2(2);
Alpha   = Alpha_n/(Gamma*Alpha_1);
Beta    = -(Beta_n + Gamma.*Beta_1 + Delta)/(Gamma.*Alpha_1);

InfoX.X2toX1 = [Alpha, Beta];
InfoX.X1toX2 = [1./Alpha, -Beta./Alpha];


%X12 = InfoX.X1toX2(1).*X1 + InfoX.X1toX2(2);
%X21 = InfoX.X2toX1(1).*X2 + InfoX.X2toX1(2);
%S12 = interp1(X1,S1,X21);
%
%plot(X21,S12,'b-');
%hold on;
%plot(X2,S2,'r-');
