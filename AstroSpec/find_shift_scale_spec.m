function [Shift,Scale,Corr,Np,RMS,Res,ScaledRef]=find_shift_scale_spec(Spec,Ref,ShiftMethod,ScaleMethod,SigmaClip,Poly)
%--------------------------------------------------------------------------
% find_shift_scale_spec function                                 AstroSpec
% Description: Given a spectrum and a reference spectrum [Wavelength, Flux]
%              find the shift (a constant added to the data) and/or a
%              scaling (multiply the data by a constant) that relats the
%              data and the reference by minimizing the rms.
%              Note that the shift is applaied first.
%              Alternatively, the program can find a factor as function
%              of wavelength that such: Ref*Factor = Spec.
%              This is done by fitting a polynomial to the log of the ratio
%              between the spectra and reference.
% Input  : - Spectrum [Wavelength, Intensity]
%          - Reference spectrum [Wavelength, Intensity]
%          - Shift method:
%            'none'   
%            'mean'
%            'median'   (default).
%            'std'
%            'min'
%            'max'
%            'fit'      - find shift and scale using a LSQ wavelength independent fit.
%          - Scale method:
%            'none'   
%            'mean'
%            'median'    (default).
%            'std'
%            'range'
%            'min'
%            'max'
%            'fit'       - find shift and scale using a LSQ wavelength independent fit.
%          - Sigma clipping (for the fit option) in std units (default is NaN).
%            if NaN then no sigma clipping.
%          - Match Spec and Ref by devision and polynomial fit.
%            If 0 or not given then use Shift and Scale method,
%            else use this method instead of Shift and Scale methods. Default is 0.
% Output : - Shift, or alternativel if Poly>0, then:
%            return [Wavelength, Factor],
%            in which wavelength is taken from the Ref, and Ref*Factor = Spec.
%          - Scale
%          - Correlation between spectra and reference in overlap region.
%          - The number of spectral points in the Spectra within the overlap region.
%          - RMS of spectra and reference residual in overlap region.
%          - Vector of residuals between spectra and reference in overlap region.
%          - Scaled reference spectra
% See also: scale_spectrum.m
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

ColW         = 1;
ColI         = 2;
ColE         = 3;
InterpMethod = 'linear';

if (nargin==2),
   ShiftMethod = 'median';
   ScaleMethod = 'median';
   SigmaClip   = NaN;
   Poly        = 0;
elseif (nargin==3),
   ScaleMethod = 'median';
   SigmaClip   = NaN;
   Poly        = 0;
elseif (nargin==4),
   SigmaClip   = NaN;
   Poly        = 0;
elseif (nargin==5),
   Poly        = 0;
elseif (nargin==6),
   % do nothing
else
   error('Illegal number of input arguments');
end   

%-------------------------------
%--- Interpolate Ref to Spec ---
%-------------------------------
OrigRef = Ref;
NewRefI = interp1(Ref(:,ColW),Ref(:,ColI),Spec(:,ColW),InterpMethod);
Ref     = [Spec(:,ColW), NewRefI];

%------------------
%--- Remove NaN ---
%------------------
I       = find(isnan(Spec(:,ColI))==0 & isnan(Ref(:,ColI))==0);
Ref     = Ref(I,:);
Spec    = Spec(I,:);


if (Poly==0),
   switch ShiftMethod
    case 'none'
       Shift = 0;
    case 'mean'
       Shift = mean(Ref(:,ColI)) - mean(Spec(:,ColI));
    case 'median'
       Shift = median(Ref(:,ColI)) - median(Spec(:,ColI));
    case 'std'
       Shift = std(Ref(:,ColI)) - std(Spec(:,ColI));
    case 'min'
       Shift = min(Ref(:,ColI)) - min(Spec(:,ColI));
    case 'max'
       Shift = max(Ref(:,ColI)) - max(Spec(:,ColI));
    case 'fit'
       % do nothing
    otherwise
       error('Unknown ShiftMethod Option');
   end
   
   switch ScaleMethod
    case 'none'
       Scale = 1;
    case 'mean'
       Scale = mean(Ref(:,ColI))./mean(Shift+Spec(:,ColI));
    case 'median'
       Scale = median(Ref(:,ColI))./median(Shift+Spec(:,ColI));
    case 'std'
       Scale = std(Ref(:,ColI))./std(Shift+Spec(:,ColI));
    case 'range'
       Scale = range(Ref(:,ColI))./range(Shift+Spec(:,ColI));
    case 'min'
       Scale = min(Ref(:,ColI))./min(Shift+Spec(:,ColI));
    case 'max'
       Scale = max(Ref(:,ColI))./max(Shift+Spec(:,ColI));
    case 'fit'
       switch ShiftMethod
        case 'fit'
           SizeSpec = size(Spec,1);
           SizeRef  = size(Ref,1);
   
           Par      = [ones(SizeRef,1),Spec(:,ColI)]\Ref(:,ColI);
           Shift    = Par(1);
           Scale    = Par(2);
   
           if (isnan(SigmaClip)==1),
              % no sigma clipping
           else
   
              Res     = Spec(:,ColI) - (Shift + Scale.*Ref(:,ColI));
              StdRes  = std(Res); 
              Isc     = find(abs(Res)<=StdRes.*SigmaClip);
              Ref     = Ref(Isc,:); 
              SizeRef = size(Ref,1);
              Spec    = Spec(Isc,:); 
   
              Par     = [ones(SizeRef,1),Spec(:,ColI)]\Ref(:,ColI);
              Shift   = Par(1);
              Scale   = Par(2);
   
           end
   
        otherwise
           error('ShiftMethod & ScaleMethod must be: fit');
       end
    otherwise
       error('Unknown ScaleMethod Option');
   end

   if (nargout>2),
      ScaledRef = [Ref(:,ColW), Shift + Scale.*Ref(:,ColI)];
   end
else
   %--- Polynomial fititing of spectra devision ---
   % polynomial scaling
   % interpolate Template to Spec
   NewRefI = interp1(Ref(:,ColW),Ref(:,ColI),Spec(:,ColW),InterpMethod);
   Factor   = Spec(:,ColI)./NewRefI;
   
   % poly fit Factor
   In0            = find(Factor>0 & isnan(Factor)==0);
   Par            = polyfit(Spec(In0,ColW),log10(Factor(In0)),Poly);
   RefFactor      = 10.^(polyval(Par,OrigRef(:,ColW)));
   Shift          = [OrigRef(:,ColW), RefFactor];
   Scale          = NaN;

   if (nargout>2),
      RefFactor      = 10.^(polyval(Par,Ref(:,ColW)));
      ScaledRef = [Ref(:,ColW), Ref(:,ColI).*RefFactor];
   end
end




if (nargout>2),
   Res       = Spec(:,ColI) - ScaledRef(:,ColI);
   RMS       = std(Res);
   MeanRef   = mean(ScaledRef(:,ColI));
   MeanSpec  = mean(Spec(:,ColI));
   Np        = size(Spec,1);
   Corr      = sum((ScaledRef(:,ColI) - MeanRef).*(Spec(:,ColI) - MeanSpec))./(Np.*std(ScaledRef(:,ColI)).*std(Spec(:,ColI)));
end
