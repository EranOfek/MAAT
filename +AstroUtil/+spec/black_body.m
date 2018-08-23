function [Il,In,IlA,ImJy,Ip]=black_body(T,W,Type)
% Black body spectrum
% Package: AstroUtil.spec
% Description: Black body spectrum.
%              OBSOLETE: Use AstSpec.blackbody instead.
% Input  : - Temperature [K].
%          - Vector of wavelength [Ang].
%          - Calculation type:
%            'P'  - planck formula, default.
%            'RJ' - Rayleigh-Jeans approximation
%            'W'  - Wein spectrum.
% Output : - Emittance [erg/sec/cm^2/cm(lambda)]
%          - Emittance [erg/sec/cm^2/Hz]
%          - Emittance [erg/sec/cm^2/Ang(lambda)]
%          - Emittance [mJy] (i.e., (erg/sec/cm^2/Hz)/1e-26) 
%          - Number of photons [photons/sec/cm^2/Ang(lambda)]
% See also: AstSpec.blackbody
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jan 2003
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:  [Il,In,IlA,ImJy,Ip]=AstroUtil.spec.black_body(10000,logspace(1,4,100).');
% Reliable: 1
%--------------------------------------------------------------------------
Def.Type = 'P';
if (nargin==2)
   Type = Def.Type;
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end

h      = 6.6261e-27;      % = get_constant('h','cgs');          % Planck constant [cgs] 
c      = 29979245800;     % = get_constant('c','cgs');          % speed of light [cm]
k      = 1.380648813e-16; % = get_constant('kB','cgs');         % Boltzmann constant [cgs]

Lam    = W.*1e-8;            % convert Ang to cm
Nu     = c./Lam;             % convert wavelength to frequency [Hz]


switch lower(Type)
 case 'p'
    %----------------------
    %--- planck formula ---
    %----------------------

    %--- wavelength units ---
    % emittence [erg/sec/cm^2/cm(lambda)]
    Il = 2.*pi.*h.*c.^2.*Lam.^(-5) ./ (exp(h.*c./(Lam.*k.*T)) - 1);
    
    % Photon distribution law
    % 2.*pi.*c.*Lam.^(-4) ./ (exp(h.*c./(Lam.*k.*T)) - 1);
    
    
    %--- frequency units ---
    % emittence [erg/sec/cm^2/Hz]
    if (nargout>1)
        In = 2.*pi.*h.*Nu.^3.*c.^(-2) ./ (exp(h.*Nu./(k.*T)) - 1);
    end
    
    % Photon distribution law
    % 2.*pi.*c.^(-2).*Nu.^2 ./ (exp(h.*Nu./(k.*T)) - 1);
    
    if (nargout>2)
       IlA  = Il.*1e-8;
       ImJy = In./1e-26;   
    
       EnConv = convert.energy('A','erg',W);
       Ip = IlA./EnConv;
    end
 case 'w'
    %----------------------
    %--- planck formula ---
    %----------------------

    %--- wavelength units ---
    % emittence [erg/sec/cm^2/cm(lambda)]
    Il = 2.*pi.*h.*c.^2.*Lam.^(-5) ./ (exp(h.*c./(Lam.*k.*T)));
    
    % Photon distribution law
    % 2.*pi.*c.*Lam.^(-4) ./ (exp(h.*c./(Lam.*k.*T)));
    
    
    %--- frequency units ---
    % emittence [erg/sec/cm^2/Hz]
    if (nargout>1)
        In = 2.*pi.*h.*Nu.^3.*c.^(-2) ./ (exp(h.*Nu./(k.*T)));
    end
    
    % Photon distribution law
    % 2.*pi.*c.^(-2).*Nu.^2 ./ (exp(h.*Nu./(k.*T)));
    
    if (nargout>2)
       IlA  = Il.*1e-8;
       ImJy = In./1e-26;   
    
       EnConv = convert.energy('A','erg',W);
       Ip = IlA./EnConv;
    end

 case 'rj'
    %-----------------------------------
    %---Rayleigh-Jeans approximation ---
    %-----------------------------------

    %--- wavelength units ---
    % emittence [erg/sec/cm^2/cm(lambda)]
    Il = 2.*pi.*k.*T.*c.*Lam.^(-4);
    
    %--- frequency units ---
    % emittence [erg/sec/cm^2/Hz]
    if (nargout>1)
        In = 2.*pi.*k.*T.*(Nu./c).^2;
    end
    
    if (nargout>2)
       IlA  = Il.*1e-8;
       ImJy = In./1e-26;   
    
       EnConv = convert.energy('A','erg',W);
       Ip = IlA./EnConv;
    end


 otherwise
    error('Unknown Type option');
end
