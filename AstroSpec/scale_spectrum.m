function NewSpec=scale_spectrum(Spec,ScalePar)
%--------------------------------------------------------------------------
% scale_spectrum function                                        AstroSpec
% Description: Scale spectrum by shift and stretch or a wavelength
%              dependent factor (See also: find_shift_scale_spec.m).
% Input  : - Spectrum to scale [Wavelength, Intensity, Error].
%            Error is optional.
%          - Scale parameters:
%            if two element vector is given then: [Shift Stretch]
%            else regarded as a matrix of [wavelength, Factor]
%            in which factor is multiplied by the intensity
%            to get the new spectrum.
% Output : - Scaled spectrum [Wavelength, Intensity, Error].
%            NewSpec = (Spec+Shift)*Stretch.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: find_shift_scale_spec.m
% Reliable: 1
%--------------------------------------------------------------------------
ColW = 1;
ColI = 2;
ColE = 3;
InterpMethod = 'linear';

NewSpec = Spec;
Ncol    = size(Spec,2);

if (numel(ScalePar)==2),
   NewSpec(:,ColI) = (NewSpec(:,ColI) + ScalePar(1)).*ScalePar(2);
   if (Ncol>2),
      NewSpec(:,ColE) = (NewSpec(:,ColE) + ScalePar(1)).*ScalePar(2);
   end
else
   % interpolate stretch factor
   Factor = interp1(ScalePar(:,1),ScalePar(:,2),Spec(:,ColW),InterpMethod);
   NewSpec(:,ColI) = NewSpec(:,ColI).*Factor;
   if (Ncol>2),
      NewSpec(:,ColE) = NewSpec(:,ColE).*Factor;
   end
end
