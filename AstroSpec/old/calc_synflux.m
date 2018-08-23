function [TotFlux,Unc,Mag_AB,Mag_Vega] = calc_synflux(ObjSpectra,Z,Filter,InterpMethod);
%----------------------------------------------------------------
% calc_synflux function                                AstroSpec
% Description: Calculate synthetic flux in a given filter from
%              a spectrum. The results are non zero-point calibrated.
% Input  : - Object spectra, two column matrix
%            [Wavelength(\AA), flux].
%          - Applay redshift to spectra, from rest frame to observed...
%          - Filter name or filter data matrix.
%            If two column matrix is given, then taken
%            as [wavelength(\AA), Transmission].
%            else give filter name (see get_filter_old.m for all optoins): 
%            {'hH' | 'hV' | 'hR' | 'hI' |
%            {'jU' | 'jB' | 'jV' | 'jR' | 'jI' |
%             'jK' | 'Su' | 'Sg' | 'Sr' | 'Si' | 'Sip' | 'Sz'}
%          - Interpolation method:
%            'linear'  - linear interpolation (default).
%            'nearest' - nearest neighbor interpolation
%            'spline'  - cubic spline interpolation
%            'cubic'   - cubic interpolation
% Output : - Total flux.
%          - Uncertainty due to non-overlap between the filter
%            and spectrum, in units of flux within filter.
%            Values near 0 indicating almost complete overlap.
% Remarks: Return un-calibrated magnitude!
%          Use zero_point to calibrate magnitude.
% Tested : Matlab 5.3
%     By : Eran O. Ofek             August 2000   
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%----------------------------------------------------------------
if (nargin==3),
   InterpMethod = 'linear';
elseif (nargin==4),
   % no default.
else
   error('Illegal number of input argumemts');
end

if (isstr(Filter)==1),
   [FilterData, ZP_AB, ZP_Vega] = get_filter_old(Filter);
else
   FilterData = Filter;
end


NewSpec = shift_spec(ObjSpectra,Z,'f','wi');

[NewList1,NewList2]=eq_sampling(NewSpec,FilterData,[],InterpMethod);

% Uncertainty
Unc = 1-trapz(NewList2(:,1),NewList2(:,2))./trapz(FilterData(:,1),FilterData(:,2));

TotFlux  = trapz(NewList1(:,1), NewList1(:,2).*NewList2(:,2));

Mag_AB   = -2.5.*log10(TotFlux) + ZP_AB;
Mag_Vega = -2.5.*log10(TotFlux) + ZP_Vega;
