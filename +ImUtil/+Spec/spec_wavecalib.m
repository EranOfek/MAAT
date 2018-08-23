function spec_wavecalib(Spec,RefSpec,Method,varargin);
%-----------------------------------------------------------------------------
% spec_wavecalib function                                              ImSpec
% Description: Given a two column matrix representing a sky or lamp spectrum
%              attempt to wavelength calibrate the spectrum.
% Input  : - Two column matrix representing a sky or lamp spectrum,
%            [pixel, intensity].
%          - calibrated reference spectrum by which to calibrate the
%            input spectrum.
%            If empty matrix, then calibrate the spectrum manualy.
%            (i.e., the user will be prompted to select line in the
%            spectrum and specify their wavelength.
%            Alternatively, this can be a two column matrix containing
%            a calibrated reference spectrum [wavelength[A], intensity].
%            Or this can be a string specifying the refernce wavelength
%            calibration spectrum name:
%            'all'  - attempt to match all reference spectra in library.
%                     see add_wavecalib_spec.m to add reference spectrum
%                     to library (default).
%            or any other spectra name (see get_wavecalib_spec.m).
%          - Method by which to perform the calibration:
%            'xcorr' - cross correlate spectra (default).
%            'man'   - manual calibration (select lines manualy).
%          * Arbitrary number of pairs of input arguments:
%            ...,'keyword','value',...
%            Optional keywords are:
%
% Output : - 
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    August 2010
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------------------

Def.RefSpec = 'all';
Def.Method  = 'xcorr';
if (nargin==1),
   RefSpec  = Def.RefSPec;
   Method   = Def.Method;
elseif (nargin==2),
   Method   = Def.Method;
else
   % do nothing
end

DefV.

InPar = set_varargin_keyval(DefV,'y','use',varargin{:});



