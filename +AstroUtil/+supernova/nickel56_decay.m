function [E,E_Ni,E_Co]=nickel56_decay(Time)
% Energy production of Nickel56->Cobalt->Iron
% Package: AstroUtil.supernova
% Description: Calculate the energy production of Nicel56->Cobalt->Iron
%              radioactive decay as a function of time.
% Input  : - Time [day].
% Output : - Total energy production [erg g^-1 s^-1] for each Time.
%          - Energy production from the Nickel decay
%            [erg g^-1 s^-1] for each Time.
%          - Energy production from the Cobalt decay
%            [erg g^-1 s^-1] for each Time.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Aug 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference : Kulkarni (2006), Eq. 14, 44.
% Example: [E,E_Ni,E_Co]=AstroUtil.supernova.nickel56_decay([1;2]);
% Reliable: 2
%--------------------------------------------------------------------------

L_Ni = 1./8.8;
L_Co = 1./111.09;
E_Ni = 3.9e10.*exp(-L_Ni.*Time);
E_Co = 7e9.*(exp(-L_Co.*Time) - exp(-L_Ni.*Time));

E = E_Ni + E_Co;
