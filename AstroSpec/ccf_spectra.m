function [CC,Npt,RedshiftVec]=ccf_spectra(Spectra,Template,RedshiftVec,ShiftMethod,ScaleMethod,SigmaClip,Poly)
%--------------------------------------------------------------------------
% ccf_spectra function                                           AstroSpec
% Description: Cross-Correlate two spectra. For each redshift, the
%              correlation is calculated by applying a redshift to the
%              second spectrum (Template), and correlate it with the
%              first spectrum.
% Input  : - Specta [Wavelength, specific flux].
%          - second spectra (Template) [Wavelength, specific flux].
%          - Redshift vector (default is: [0]).
%          - Match spectra before coorelation, Shift method
%            {'none' | 'mean' | 'median' | 'std' | 'min' | 'max' | 'fit' | 'old'},
%            default is 'none' - see find_shift_scale_spec.m for details.
%            Note that the matching is done only in the overlaping region
%            unless the 'old' option is specified in which the old 
%            algorithm using nanmean allover the spectrum.
%          - Match spectra before coorelation, Scale method
%            {'none' | 'mean' | 'median' | 'std' | 'min' | 'max' | 'fit'},
%            default is 'median' - see find_shift_scale_spec.m for details.
%          - Sigma clipping (for the fit option) in std units (default is NaN).
%            if NaN then no sigma clipping.
%          - Match Spec and Ref by devision and polynomial fit.
%            If 0 or not given then use Shift and Scale method,
%            else use this method instead of Shift and Scale methods. 
%            Default is 0. see find_shift_scale_spec.m for details.
% Output : - Correlation coef. vector for each redshift.
%          - The number of points in the overlap region.
%          - Redshift vector.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 1
%--------------------------------------------------------------------------
SampRate = 1;
Method   = 'linear';
if (nargin==2),
   RedshiftVec = [0];
   ShiftMethod = 'none';
   ScaleMethod = 'median';
   SigmaClip   = NaN;
   Poly        = 0;
elseif (nargin==3),
   ShiftMethod = 'none';
   ScaleMethod = 'median';
   SigmaClip   = NaN;
   Poly        = 0;
elseif (nargin==4),
   ScaleMethod = 'median';
   SigmaClip   = NaN;
   Poly        = 0;
elseif (nargin==5),
   SigmaClip   = NaN;
   Poly        = 0;
elseif (nargin==6),
   Poly        = 0;
elseif (nargin==7),
   % do nothing
else
   error('Illegal number of input arguments');
end

N     = length(RedshiftVec);
CC    = zeros(size(RedshiftVec));
Npt   = zeros(size(RedshiftVec));

for I=1:1:N,
   Z       = RedshiftVec(I);
   RedSpec = shift_spec(Template, Z, 'f','wi');

   switch ShiftMethod
    case 'old'
       %--- old correlation ---
       [NewList1,NewList2]=eq_sampling(Spectra, RedSpec, SampRate, Method);

        Np       = length(NewList1(:,1));
        CC(I)    = nansum((NewList1(:,2) - nanmean(NewList1(:,2))).*  ...,
                          (NewList2(:,2) - nanmean(NewList2(:,2))))./ ...,
                          (Np.*nanstd(NewList1(:,2)).*nanstd(NewList2(:,2)));
        Npt(I)   = nanmax(NewList1(:,1)) - nanmin(NewList1(:,1));

    otherwise
       %--- new correlation ---
       [Shift,Scale,Corr,Np,RMS,Res,ScaledRef] = find_shift_scale_spec(Spectra,RedSpec,ShiftMethod,ScaleMethod,SigmaClip,Poly);

        CC(I)              = Corr;
        Npt(I)             = Np; 

   end
end



