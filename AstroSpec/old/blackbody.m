function [AstS]=blackbody(VecT,VecW,UnitsOut,UnitsWave)
%--------------------------------------------------------------------------
% blackbody function                                             AstroSpec
% Description: Calculate black body (planck) spectrum.
% Input  : - Array of temperatures [K]. Calculate a planck spectrum for
%            each temperature.
%          - Vector of wavelength. Default is (1000:100:25000)';
%          - Output units. See convert.flux for options.
%            Default is 'cgs/A'.
%          - Input wavelength units. See convert.units for options.
%            Default is 'ang'.
% Output : - AstSpec array with black body spectrum.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstS=blackbody([5000;10000]);
% Reliable: 2
%--------------------------------------------------------------------------

Def.VecW      = (1000:100:25000).';
Def.UnitsOut  = 'cgs/A';
Def.UnitsWave = 'Ang';
if (nargin==1),
    VecW      = Def.VecW;
    UnitsOut  = Def.UnitsOut;
    UnitsWave = Def.UnitsWave;
elseif (nargin==2),
    UnitsOut  = Def.UnitsOut;
    UnitsWave = Def.UnitsWave;
elseif (nargin==3),
    UnitsWave = Def.UnitsWave;
elseif (nargin==4),
    % do nothing
else
    error('Illegal number of input arguments');
end

h      = 6.6261e-27;      % = get_constant('h','cgs');          % Planck constant [cgs] 
c      = 29979245800;     % = get_constant('c','cgs');          % speed of light [cm]
k      = 1.380648813e-16; % = get_constant('kB','cgs');         % Boltzmann constant [cgs]

Lam = VecW.*convert.units(UnitsWave,'cm');   % convert wavelength to cm
%Lam    = W.*1e-8;            % convert Ang to cm
Nu     = c./Lam;             % convert wavelength to frequency [Hz]

Nt = numel(VecT);
AstS = astspecdef(size(VecT));
for It=1:1:Nt,
    % for each temperature
    AstS(It).Wave = VecW;
    
    %----------------------
    %--- planck formula ---
    %----------------------
    % emittence [erg/sec/cm^2/Ang(lambda)]
    Flambda = 1e-8.*2.*pi.*h.*c.^2.*Lam.^(-5) ./ (exp(h.*c./(Lam.*k.*VecT(It))) - 1);
 
    AstS(It).Int = convert.flux(Flambda,'cgs/A',UnitsOut,Lam,'cm');
    AstS(It).WaveUnits = UnitsWave;
    AstS(It).IntUnits  = UnitsOut;
    AstS(It).source    = 'blackbody.m';
    AstS(It).ObjName   = sprintf('Planck spectrum T=%f',VecT(It));
    AstS(It).z         = 0;
end
