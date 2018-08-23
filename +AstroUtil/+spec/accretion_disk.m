function [Spec,VecR,VecT,Rs]=accretion_disk(M,Mdot,Rin,Rout,Wave,Nstep)
% theoretical spectrum of a optically thick, thin accretion disk
% Package: AstroUtil.spec
% Description: Calculate the theoretical spectrum of a optically thick,
%              thin accretion disk given the mass of
%              the central object, disk radius, and the accretion rate.
% Input  : - Mass of central object [Solar Mass].
%          - Accreation rate (Mdot) [Solar Mass per Year].
%          - Inner cutof radius of accretion disk [R_{s}],
%            default is 10.
%            If empty the use default.
%          - Outer cutof radius of accretion disk [R_{s}],
%            default is 1000.
%            If empty the use default.
%          - Wavelength in which to calculate Flux [Ang], default is
%            logspace(0,4,1000)
%            If empty the use default.
%          - Number of steps in accretion disk radius, default is 1000.
%            If empty the use default.
% Output : - Total emittance [erg/sec/cm^2/cm(lambda)] from accretion
%            disk as function of wavelength, [Wavelengh[Ang], Lum].
%          - Vector of radii in the accretion disk [cm].
%          - Vector of temperature corresponding to eacg radius [K].
%          - The R_s radius of the central mass.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Spec,VecR,VecT,Rs]=AstroUtil.spec.accretion_disk(1e6,0.1,10,100);
% Reliable: 2
%--------------------------------------------------------------------------
Rin_Def    = 10;
Rout_Def   = 1000;
Wave_Def   = logspace(0,4,1000).';
Nstep_Def  = 1000;

if (nargin==2)
   Rin     = Rin_Def;
   Rout    = Rout_Def;
   Wave    = Wave_Def;
   Nstep   = Nstep_Def;
elseif (nargin==3)
   Rout    = Rout_Def;
   Wave    = Wave_Def;
   Nstep   = Nstep_Def;
elseif (nargin==4)
   Wave    = Wave_Def;
   Nstep   = Nstep_Def;
elseif (nargin==5)
   Nstep   = Nstep_Def;
elseif (nargin==6)
   % do nothing
else
   error('Illegal number of input arguments');
end
if (isempty(Rin))
   Rin = Rin_Def;
end
if (isempty(Rout))
   Rout = Rout_Def;
end
if (isempty(Wave))
   Wave = Wave_Def;
end
if (isempty(Nstep))
   Nstep = Nstep_Def;
end


%--- Constants ---
Msun     = constant.SunM;      % solar mass [gram]
SecYear  = 365.25.*86400; %get_constant('JY','cgs');        % sec in year
G        = constant.G;         % gravitational constant [cgs]
Sig      = constant.sigma;     % Stefan-Boltzmann constant [cgs]
c        = constant.c;         % speed of light [cm]

M    = M.*Msun;                   % convert solar mass to grams
Mdot = Mdot.*Msun./SecYear;       % convert solar/year to gram/sec

Rs   = 2.*G.*M./(c.^2);
Rin  = Rs.*Rin;
Rout = Rs.*Rout;

Sum    = zeros(size(Wave));

% integrated spectrum
StepR = (Rout - Rin)./Nstep;
VecR = (Rin:StepR:(Rout-StepR)).';
I = 0;
VecT = zeros(size(VecR));
for R=Rin:StepR:(Rout-StepR)
   R       = R + 0.5.*StepR;
   Area    = 2.*pi.*R.*StepR;

   % temperature (T) as function of radius (R)
   T       = (3.*G.*M.*Mdot./(8.*pi.*Sig.*R.^3)  .* (1 - sqrt(Rin./R)) ).^0.25;
   I = I + 1;
   VecT(I) = T;
   [Il,~] = AstroUtil.spec.black_body(T,Wave);

   Sum     = Sum + Il.*Area;
end

Spec = [Wave, Sum];


