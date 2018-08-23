function T=brightness_temp(S,Freq,Dist,R)
%--------------------------------------------------------------------------
% brightness_temp function                                       AstroSpec
% Description: Calculate the brightness temperature.
% Input  : - Specific flux [erg cm^-2 s^-1 Hz^-1]
%            (note that 1mJy = 1e-26 ergs cm^-2 s^-1 Hz^-1).
%          - Frequency [Hz].
%          - Luminosity distance [pc].
%          - Source radius [cm].
% Output : - Brightness temperature [K].
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Nov 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: T=brightness_temp(1e-26,1.4e9,10,1e6);
% Reliable: 1
%--------------------------------------------------------------------------

k  = get_constant('kB','cgs');
c  = get_constant('c','cgs');
pc = get_constant('pc','cgs');
Dist = Dist.*pc;

T = Dist.^2.*0.5.*S.*(c./Freq).^2./(pi.*k.*(R.^2));

