function Res=sedm_corr_waveflexure(DY,HA,Dec,PixScale)
%--------------------------------------------------------------------------
% sedm_corr_waveflexure function                                      SEDM
% Description: Estimate the flexure in the dispersion direction
%              (wavelength) in the SEDM.
% Input  : - Vector or scalar of slope Y-position offset [pix] relative to
%            segmentation map taken at zenith.
%          - Hour Angle of observations [deg].
%          - Declination of observations [deg].
%          - Approximate pixel scale [Ang/Pix]. Default is 25.
% Output : - Approximate correction due to flexure in wavelength.
%            This is measured relative to an arc that was taken in the
%            zenith direction.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Aug 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=sedm_corr_waveflexure(0,0,0);
% Reliable: 
%--------------------------------------------------------------------------


Def.PixScale = 25;  % [Ang/pix]
if (nargin==3),
    PixScale = Def.PixScale;
end

% 
Res = -9.76170735e-2 ...
      -5.40857690e-1.*DY ...
      -1.42770036e-2.*HA ...
      -8.24304496e-3.*Dec ...
      +9.89255017e-3.*DY.^2 ...
      +3.53396833e-4.*DY.*HA ...
      +2.31293635e-3.*DY.*Dec ...
      +7.37365997e-5.*HA.^2 ...
      +8.83937326e-5.*HA.*Dec ...
      -6.34635265e-6.*Dec.^2;
  
  Res = Res.*PixScale;