function IllumMat=spec_illum(Mat,SemiBinX,FitMethod,Order,CollapseMethod);
%--------------------------------------------------------------------------
% spec_illum function                                               ImSpec
% Description: Create an illumination calibration image for
%                         long slit spectrum.
%                         The illumination calibration is performed by
%                         fitting polynomial functions across
%                         the slit (the slit  profiles) at conscutive
%                         dispersion points (with binning) and
%                         normalizing each fitted function to unity at
%                         the center of the slit.
% Input  : - Image matrix, x-axis is the dispersion axis.
%          - Half width of x-axis (dispersion-axis) bin size.
%          - Fit method:
%            'Poly' - polynomial fit
%            'runmean' - running mean with a gaussian window function
%          - Order of fit method.
%            In case of polynomail fit then this is the order of the polynomial.
%            In case of runmean this is the width of the gaussian.
%          - avergaing method in each bin size
%            (averaging in dispersion direction):
%            {'median' | 'mean'}, default is 'median'.
% Output : - Illumination matrix.
% Tested : Matlab 6.5
%     By : Eran O. Ofek                 September 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------------

CollapseMethodDef = 'median';

if (nargin==4),
   CollapseMethod = CollapseMethodDef;
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end

[Ny,Nx] = size(Mat);
VecY    = [1:1:Ny].';
VecX    = [1:1:Nx].';
CenterY = floor(0.5.*Ny);

% inilialize IllumMat (Illumination image)
IllumMat = zeros(Ny,Nx);

for I=1:1:Nx,
   CurrentX = VecX(I);
   % lower boundry of SubX region
   MinSubX  = max([1; CurrentX-SemiBinX]);
   % upper boundry of SubX region
   MaxSubX  = min([Nx; CurrentX+SemiBinX]);
   % SubX region
   SubVecX  = [MinSubX:1:MaxSubX].';

   SubMat   = Mat(VecY,SubVecX);

   %--- collapse SubMat over dispersion axis ---
   switch CollapseMethod
    case 'median'
       CollapseSubMat = median(SubMat,2);
    case 'mean'
       CollapseSubMat = mean(SubMat,2);
    otherwise
       error('Unknwon CollpaseMethod option');
   end

   switch FitMethod
    case 'Poly'
       %--- fit polynomial to CollapseSubMat ---
       PolyPar = polyfit(VecY,CollapseSubMat,Order);
       SmoothCollapse  = polyval(PolyPar,VecY);
    case 'runmean'
       %--- running mean ---
       RunMean = runmean([VecY,CollapseSubMat],'g',Order,'f');
       SmoothCollapse  = RunMean(:,2);
    otherwise
       error('Unkown FitMethod option');
   end

   IllumMat(VecY,CurrentX) = SmoothCollapse./SmoothCollapse(CenterY);
end

   
