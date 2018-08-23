function [SkySubSpec,SkyVal,SkyRMS]=spec_skysub1(ExtractedSpec,BackReg,varargin);
%---------------------------------------------------------------------------
% spec__skysub function                                              ImSpec
% Description: Subtract the sky from a 2-dimensional spectrum.
% Input  : - Matrix containing a 2-dimensional spectrum from which to
%            subtract the sky. In the default mode, the dispersion axis
%            is along the x-axis.
%          - Background region in which to fit the background. This is
%            a 2x2 matrix in which the first column indicate the lower
%            and upper position of the left background region, and
%            the second column indicate the lower and upper position of
%            the right background region.
%            If empty matrix, or not given then use interactive mode.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...,
%            where keywords can be one of the followings:
%            'SkyMet'   - Sky fitting method (default is  'poly2').
%                         'mean'     - mean of sky region.
%                         'median'   - median of sky region.
%                         'poly#'    - polynomial fit with optional CR
%                                      cleaning of sky region. The degree
%                                      of polynomial is given in the # sign
%                                      (e.g., poly2').
%            'Nrej'     - If the sky fitting algorith is using a cleaning
%                         scheme, then set [Nlow, Nhigh], which are the
%                         number of lowest/highest N pixels to reject
%                         from the fit after the first iteration.
%                         Default is [1 3].
%            'ColMet'   - Collapse method {'mean' | 'median'},
%                         default is 'median'.
%            'DispAxis' - Dispersion axis {'x'|'y'}, default is 'x'.
% Output : - Sky subtracted spectrum, in which the disperion axis is
%            along the x-axis (regardless of the imput).
%          - Vector of best fit sky level for each position in the input
%            spectrum. The disperion axis is along the x-axis (regardless
%            of the input).
%          - Vector of sky background RMS based on the best fitted sky.
%            The disperion axis is along the x-axis (regardless of the
%            input).
% Tested : Matlab 5.3
%     By : Eran O. Ofek                     April 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%---------------------------------------------------------------------------
DefBackReg        = [];
if (nargin==1),
   BackReg        = DefBackReg;
else
   % do nothing
end


% set default parameters
SkyMet      = 'poly2';
Nrej        = [1 3];
ColMet      = 'median';
DispAxis    = 'x';

Narg  = length(varargin);
for Iarg=1:2:Narg-1,
   switch lower(varargin{Iarg})
    case 'skymet'
       SkyMet     = varargin{Iarg+1};
    case 'nrej'
       Nrej       = varargin{Iarg+1};
    case 'colmet'
       ColMet     = varargin{Iarg+1};
    case 'dispaxis'
       DispAxis   = varargin{Iarg+1};
    otherwise
       error('Unknown keyword option');
   end
end

Nlow    = Nrej(1);
Nhigh   = Nrej(2);

switch lower(DispAxis)
 case 'x'
    % do nothing
 case 'y'
    % rotate
    ExtractedSpec    = ExtractedSpec.';
 otherwise
    error('Unknown DispAxis option');
end
DispAxisN = 2;

[Nspat,Ndisp] = size(ExtractedSpec);
VecDisp  = [1:1:Ndisp].';
VecSpat  = [1:1:Nspat].';


switch ColMet
 case 'mean'
    Collapse = mean(ExtractedSpec,DispAxisN);
 case 'median'
    Collapse = median(ExtractedSpec,DispAxisN);
 otherwise
    error('Unknown ColMet option');
end

if (isempty(BackReg)==1),
   %------------------------
   %--- Interactive mode ---
   %------------------------
   stairs(VecSpat,Collapse,'b');
   zoom on;
   Ans = input('Zoom is on - Click any key to continue','s');

   disp(sprintf('Select left background region (two clicks)'))
   [Xl,Yl] = ginput(2);
   disp(sprintf('Select right background region (two clicks)'))
   [Xr,Yr] = ginput(2);

   BackReg = [Xl, Xr];
else
   %----------------------------
   %--- Non interactive mode ---
   %----------------------------
   Xl  = BackReg(:,1);
   Xr  = BackReg(:,2);
end



% indices of pixels in region to fit
Ireg = find((VecSpat>=min(Xl) & VecSpat<=max(Xl)) | (VecSpat>=min(Xr) & VecSpat<=max(Xr)));


switch lower(SkyMet)
 case 'mean'
    Sky    = mean(ExtractedSpec(Ireg,:),1);
    RMS    = std(ExtractedSpec(Ireg,:),0,1);
    SkyVal = ones(Nspat,1)*Sky;
    SkyRMS = ones(Nspat,1)*RMS;

 case 'median'
    Sky    = median(ExtractedSpec(Ireg,:),1);
    RMS    = std(ExtractedSpec(Ireg,:),0,1);
    SkyVal = ones(Nspat,1)*Sky;
    SkyRMS = ones(Nspat,1)*RMS;
 otherwise
    if (strcmpi(SkyMet(1:4),'poly')==1),
       % polynomial fit
       PolyDeg   = str2num(SkyMet(5:end));

       % run for each dispersion line
       SkyVal    = zeros(Nspat,Ndisp);
       SkyRMS    = zeros(Nspat,Ndisp);
       for Idisp=1:1:Ndisp,
          Par   = polyfit(VecSpat(Ireg),ExtractedSpec(Ireg,Idisp),PolyDeg);
          Sky   = polyval(Par,VecSpat);    % best fit sky value
          Diff  = ExtractedSpec(Ireg,Idisp) - Sky(Ireg);
          RMS   = std(Diff);
          if (sum(Nrej)>0),
             % bad pixel rejection is required
             [SortedDiff,SortedInd] = sort(Diff);
             Igood = SortedInd(Nlow+1:end-Nhigh);
             IregN = Ireg(Igood);   % new pixels to fit

             Par   = polyfit(VecSpat(IregN),ExtractedSpec(IregN,Idisp),PolyDeg);
             Sky   = polyval(Par,VecSpat);    % best fit sky value
             Diff  = ExtractedSpec(IregN,Idisp) - Sky(IregN);
             RMS   = std(Diff);
          end
          SkyVal(:,Idisp) = Sky;
          SkyRMS(:,Idisp) = RMS;
       end
    else
       error('Unknown SkyMet option');
    end
end

SkySubSpec = ExtractedSpec - SkyVal;
