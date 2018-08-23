function match1d(List,RefList,Tol,Deg,SigClip);
%---------------------------------------------------------------------------
% match1d function                                                   ImSpec
% Description: Given a list of 1-dimensional position (e.g., identified
%              spectral line peaks) in an arbitrary coordinate system
%              (e.g., image coordinates), and a reference list of
%              1-dimensional position in the reference coordinate system
%              (e.g., calibrated wavelength), find the transformation
%              required in order to convert the arbitrary coordinate
%              system (of the first list) to the reference coordinate
%              system. The program assumes that the transformation is 
%              nearly linear (but not exactly).
% Input  : - Vector of coordinates.
%          - Vector of reference coordinates.
%          - Tolerance for matching the coordinates, in units...
%          - Degree of polynomial transformation, default is 2.
%          - [Low, High] sigma clipping (in StD units) for fitting the
%            polynomial transformation. Default is [-Inf Inf].
% Output : -
% Tested : Matlab 7.0
%     By : Eran O. Ofek                          April 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------
DefTol     = 1
DefDeg     = 2;
DefSigClip = [-Inf Inf];

Nlf = 6;       % Number of lines to fit

if (nargin==2),
   Tol     = DefTol;
   Deg     = DefDeg;
   SigClip = DefSigClip;
elseif (nargin==3),
   Deg     = DefDeg;
   SigClip = DefSigClip;
elseif (nargin==4),
   SigClip = DefSigClip;
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(Tol)==1),
   Tol     = DefTol;
end
if (isempty(Deg)==1),
   Deg     = DefDeg;
end
if (isempty(SigClip)==1),
   SigClip = DefSigClip;
end



% prep inverted index table
MaxRange = 15000;   % Ang
Nrl = size(RefList,1);
Nsim = 10000;
K    = 0;
for I=1:1:Nsim,
   RandI  = randperm(Nrl);
   RandI  = RandI(1:Nlf);
   Pos    = sort(RefList(RandI,1));
   Range  = range(Pos);
   if (Range<MaxRange),
      K = K + 1;
      IndTable.RelPos(K,1:Nlf) = (Pos.' - Pos(1))./Range;
      IndTable.Pos(K,1:Nlf)    = [Pos.'];
      IndTable.Range(K,1)      = Range;
   end
end

N_RelPos = K;

ThreshRelPos  = 1./100;

%PeakRegion   = 2;
%LocalBackBin = 100;
%ThreshStD    = 1.5;
%ExactPeaks = find_maxs(List,PeakRegion,LocalBackBin);
%I=find(ExactPeaks(:,5)>ThreshStD);
%ExatctPeaks = ExactPeaks(I,:);

ExactPeaks = List(:,1);
Nl = size(ExactPeaks,1);
Found.Scale = [];
Nsim = 10000;
for Isim=1:1:Nsim,
   RandI = randperm(Nl);
   RandI = RandI(1:Nlf);

   Pos     = sort(ExactPeaks(RandI,1)).';
   Range   = range(Pos);
   RelPos  = (Pos - Pos(1))./Range;

   MaxDiff = max(abs(IndTable.RelPos - ones(N_RelPos,1)*RelPos),[],2);
   II      = find(MaxDiff<ThreshRelPos);

   if (isempty(II)==0),
      Found.Scale = [Found.Scale; Range./IndTable.Range(II)];
   end
end

%   IndValMin = bin_sear(IndTable.RelPos1,RelPos1-ThreshRelPos);
%   IndValMax = bin_sear(IndTable.RelPos1,RelPos1+ThreshRelPos);
%   II = [IndValMin:IndValMax].';
%   Found.Scale = [Found.Scale; Range./IndTable.Range(II)];
   % if found - keep
   % if three (not the same) found - stop
%end


[x,n]=realhist(log10(Found.Scale),[0 3 300]);
bar(x,n)
sum(n)


Nlist    = length(List);
DiffList = diff(List);

%for Ilist=1:1:Nlist,
   % select all matches
   % between the diff of the first list and
   % reference list



%Nst = length(StretchVec);
%Nsh = length(ShiftVec);
%for Ist=1:1:Nst,
%   for Ish=1:1:Nsh,

%	       List.*StretchVec(Ist) + ShiftVec(I
