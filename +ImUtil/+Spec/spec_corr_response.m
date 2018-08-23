function [CorrSpec,InvTran]=spec_corr_response(Spec,InvTran,varargin)
%--------------------------------------------------------------------------
% spec_corr_response function                                       ImSpec
% Description: Given an observed individual spectrum and the inverse of 
%              the telescope/instrument/atmospheric response return the
%              corrected spectrum outside the atmosphere in physical
%              units.
% Input  : - Observed spectrum to correct [Wave(Ang), Intensity].
%          - Either a two column matrix [[Wave(Ang), 1./Response],
%            or a structure (e.g., returned by spec_response.m) that
%            contains fields name .VecWave and .InvTranS.
%            This response describes the wavelength dependent factor 
%            which one need to multiply the spectrum corrected for
%            air mass of zero and 1 s exposure (i.e., outside the
%            atmosphere) in order to get the spectrum in physical
%            units (e.g., erg/cm^2/s/A).
%          * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'AM'    - Airmass at which the spectrum was obtained.
%                      Default is []. If empty then use the Altitude
%                      in order to calculate the airmass.
%            'Alt'   - Altitude (radians) at which the spectrum was
%                      pbtained. Default is pi./2.
%                      If both Altitude and Airmass are provided then
%                      the Airmass will be used.
%            'ExpTime'- Exposure time (s) of spectrum observations.
%                      Default is 1.
%            'Ext'   - Atmospheric extinction type. See atmospheric_ext.m
%                      for options. Default is
%                      'KPNO_atmospheric_extinction.dat'.
%                      If empty then will not apply atmospheric extinction
%                      correction.
%            'InterpMethod' - Interpolation method. See interp1.m for
%                      options. Default is 'linear'.
% Output : - Corrected spectrum in physical units [Wavelength, Flux].
%          - Actual (resampled) inverse transmission curve used.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


TranField             = 'InvTranS';
WaveField             = 'VecWave';

DefV.AM               = [];
DefV.Alt              = pi./2;
DefV.ExpTime          = 1;
DefV.Ext              = 'KPNO_atmospheric_extinction.dat';
DefV.InterpMethod     = 'linear';

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

% set AirMass
if (isempty(InPar.AM)),
    InPar.AM = hardie(pi./2-InPar.Alt);
end

% correct spectrum for exposure time
% spectrum in 1s exposure
Spec(:,2) = Spec(:,2)./InPar.ExpTime;  % ADU/s

%--- Atmospheric extinction ---
% If required, correct the observed spectrum from
% extincted to airles spectrum
if (isempty(InPar.Ext)),
    % no atmospheric extinction correction - already airless
    AirlessSpec = Spec;
else
    AirlessSpec = atmospheric_ext(Spec,InPar.AM,InPar.Ext,InPar.InterpMethod);
end

if (~isstruct(InvTran)),
    Tran.(WaveField) = InvTran(:,1);
    Tran.(TranField) = InvTran(:,2);
else
    Tran = InvTran;
end

%--- vector of wavelength ---
%VecWave = AirlessSpec(:,1);
%Nwave   = numel(VecWave);

%--- sample transmission like spectrum ---
InvTran = interp1(Tran.(WaveField),Tran.(TranField),AirlessSpec(:,1),InPar.InterpMethod,'extrap');

%--- Correct the spectrum ---
CorrSpec = AirlessSpec;
CorrSpec(:,2) = AirlessSpec(:,2).*InvTran;
