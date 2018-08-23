function [Const,SmJy,F_lambda,F_nu]=mag_jy_conv(Filter,EffW,Alpha,System,Mag);
%---------------------------------------------------------------------
% mag_jy_conv function                                      AstroSpec
% Description: Calculate conversion factor between magnitude, flux and
%              specific flux units (e.g., mJy).
%              Allows for user supplied filter and spectral slope.
% Input  : - Filter [Wave(Ang), Trans]
%            If 'Filter' contains only one line, [wavelength, width],
%            it is interpreted as a box shape filter centered on
%            wavelength with the specified width.
%            If contains a cell array then it is regarded as
%            {Filter_Group, Filter_Name}, e.g., {'SDSS','r'}.
%          - Effective wavelength [Ang] for which to calculate
%            specific flux.
%            If empty and first input is a cell array then attempt
%            to obtain the effective wavelength using get_filter.m
%          - Spectral slope (Nu^Alpha), default is 0.
%          - Magnitude system {'Vega' | 'AB'}, default is 'Vega'.
%          - Optional magnitude (Mag) to convert to mJy.
% Output : - Mag of a 1mJy source (at the effective wavelength).
%          - The specific flux in mJy [10^-26 erg cm^-2 s^-1 Hz^-1]
%            of the optional magnitude (Mag).
%            This property is integrated over the filter width.
%          - The specific flux in [erg cm^-2 s^-1 Ang^-1]
%            of the optional magnitude (Mag).
%            This property is integrated over the filter width.
%          - The specific flux in [erg cm^-2 s^-1 Hz^-1]
%            of the optional magnitude (Mag).
%            This property is integrated over the filter width.
% Tested : Matlab 7.0
%     By : Eran O. Ofek        April 2006
%    URL : http://wise-obs,tau.ac.il/~eran/matlab.html
% Example: G=get_filter('SDSS','g');
%          [Const,SmJy,F_lambda,F_nu]=mag_jy_conv(G.nT{1},4000,0,'AB',14);
%---------------------------------------------------------------------
if (nargin==2),
   Alpha  = 0;
   System = 'Vega';
elseif (nargin==3),
   System = 'Vega';
elseif (nargin==4),
   % do nothing
elseif (nargin==5),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (iscell(Filter)),
   FilterStruct = get_filter(Filter{1},Filter{2});
   Filter = FilterStruct.nT{1};

   if (isempty(EffW)),
      EffW = FilterStruct.eff_wl{1};
   end

else
   if (size(Filter,1)==1),
      % constract a box shape filter around wavelength
      % with of the filter is taken from the second column.
      Filter = [Filter(1)-0.51.*Filter(2) 0;   Filter(1)-0.5.*Filter(2)  1;  Filter(1)+0.5.*Filter(2) 1;   Filter(1)+0.51.*Filter(2) 0];
   end
end

C  = get_constant('c','cgs');
W  = [100:100:50000].';
Nu    = C./(W.*1e-8);
EffNu = C./(EffW.*1e-8);
Flux  = 1e-26.*1e-8.* (Nu./EffNu).^Alpha .* C./((W.*1e-8).^2);
Const = basic_synthetic_photometry([W,Flux],Filter,System,[],[0 1]);

if (nargin==5 & nargout>1),
   % Convert optional magnitude to specific flux
   SmJy = 10.^(0.4.*(Const-Mag));

   F_nu = SmJy.*1e-26;    % [erg cm^-2 s^-1 Hz^-1]
   C    = get_constant('c','cgs');
   Nu   = C./(EffW.*1e-8);

   F_lambda = F_nu.*Nu.^2./C .*1e-8;   % [erg cm^-2 s^-1 Ang^-1]

   % 1 Jy = 1.51e7 photons sec^-1 m^-2 (dlambda/lambda)^-1

end


%if (nargout>4 && nargin==5),
%   % Attempt to integrate flux within filter

