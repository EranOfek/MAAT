function Spec=get_gaia_synspec(Temp,Grav,Metal,Rot,OutType)
%--------------------------------------------------------------------------
% get_gaia_synspec function                                      AstroSpec
% Description: get a synthetic stellar spectrum from the local GAIA
%              spectral library. Spectra are in the range 2500-10500A
%              and 1A resolution.
%              Assuming alpha enhanement 0, and micro-turbulence 2km/s.
% Input  : - Effective temperature [K] (in range 3500-47500 K).
%          - Gravity [log g] (in range 0 to 5).
%          - Metallicity [log solar] (in rahe -2.5 to 2.5).
%          - Rotation velocity [km/s] (in the range 0 to 500km/s).
%          - Output type: 'astspec' | 'mat'. Default is 'astspec'.
% Output : - Spectrum in flux units [wavelength[Ang], Flux].
%            Flux units [erg cm^-2 s^-1 A^-1 on star]
%            Return NaN if spectrum doesn't exist or web site is down.
% Reference: http://gaia.esa.int/spectralib/spectralib1A/SpectraLib1a.cfm
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Nov 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: wget_gaia_synspec.m
% Example: Spec=get_gaia_synspec(5000,0,0,0);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<5),
    OutType = 'astspec';
end

MSDir      = Util.files.which_dir(mfilename);
DirLocation = sprintf('%s%s..%s..%s%s%s%s%s',MSDir,filesep,filesep,filesep,'data',filesep,'GAIA_SpecTemplates',filesep);
if (Metal<0),
   MetalSign = 'M';
else
   MetalSign = 'P';
end

DefaultPars = 'K2SNWNVD01F.mat';
SpecName = sprintf('%sT%05dG%02d%s%02dV%03d%s',DirLocation,round(Temp),round(Grav.*10),MetalSign,round(abs(Metal).*10),round(Rot),DefaultPars);
W    = Util.IO.load2(sprintf('%s%s',DirLocation,'GAIA_Wave1A.mat'));
Spec = [W, load2(SpecName)];

switch lower(OutType)
    case 'astspec'
        Spec=mat2spec(AstSpec,Spec,{'Wave','Int'},{'Ang','erg*cm^-2 *s^-1*Ang^-1'});
        Spec.z = 0;
        Spec.source = 'GAIA local DB';
        Spec.ObjName = sprintf('GAIA synspec T=%f, g=%f, M=%f, R=%f',Temp,Grav,Metal,Rot);
    case 'mat'
        % do nothing
    otherwise
        error('Unknown OutType option');
end
