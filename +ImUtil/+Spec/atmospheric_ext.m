function [NewSpec,ExtFactor]=atmospheric_ext(Spec,AirMass,Ext,InterpMethod)
%------------------------------------------------------------------------------
% atmospheric_ext function                                              ImSpec
% Description: Correct a spectrum for atmospheric extinction.
% Input  : - Input spectrum [Wavelength(Ang), Flux(Arbitrary)], or
%            alternatively, two element cell array in which
%            the first contains a column vector of [Wavelength(Ang)]
%            and the second contains a matirx of fluxes, in which
%            each line correspond to each line in the wavelength vector. 
%          - Airmass in which the input spectrum was obtained.
%            If instead of correcting the spectrum for extinction you
%            would like to apply an atmospheric extinction to the spectrum,
%            then use negative airmass.
%          - Extinction type. This can be a file name or a matrix
%            containing the atmospheric extinction curve. The format
%            of the file or the matrix is: 
%            [wavelength(Ang), Extinction_per_unity_airmass(mag)].
%            Examples: 'VLT'
%                      'KPNO'
%            default is: 'VLT'.
%          - Interpolation method (see interp1.m), default is 'linear'.
% Output : - Spectrum corrected for atmospheric extinction (i.e.,
%            for airmass 0). Format is identical to this of the input
%            spectrum.
%          - Vector of extinction factor in any wavelength specified in
%            the input spectrum.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [NewSpec,ExtFactor]=atmospheric_ext([5000 1; 8000 1],1);
% Reliable: 2
%------------------------------------------------------------------------------
DefExt          = 'VLT';
DefInterpMethod = 'linear';
if (nargin==2)
   Ext          = DefExt;
   InterpMethod = DefInterpMethod;
elseif (nargin==3)
   InterpMethod = DefInterpMethod;
elseif (nargin==4)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(Ext)==1)
   Ext  = DefExt;
end

if (ischar(Ext)==1)
   %Ext = load(Ext);
   Ext = cats.spec.AtmoExtinction.(Ext);
end

% Sample Ext as same as Spec:
if (iscell(Spec)==1)
   WaveSampling = Spec{1};
else
   WaveSampling = Spec(:,1);
end
[SmSpec,SmExt] = AstroUtil.spec.eq_sampling(Spec,Ext,WaveSampling,InterpMethod);


ExtFactor = 10.^(0.4.*SmExt(:,2).*AirMass);
if (iscell(Spec)==1)
   % cell format
   NewSpec{1} = SmSpec{1};
   Ncol       = size(SmSpec{2},2);
   NewSpec{2} = [SmSpec{2}.*(ExtFactor*ones(1,Ncol))];
else
   % vector format
   NewSpec   = [SmSpec(:,1), SmSpec(:,2).*ExtFactor];
end
