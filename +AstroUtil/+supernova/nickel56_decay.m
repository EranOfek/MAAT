function [Res]=nickel56_decay(Time,t0,n)
% Energy production of Nickel56->Cobalt->Iron
% Package: AstroUtil.supernova
% Description: Calculate the energy production of Nicel56->Cobalt->Iron
%              radioactive decay as a function of time.
% Input  : - Time [day].
%          - t0. Default is 30 days.
%          - n. Default is 3.
% Output : A structure with the following fields:
%          'E' - Total energy [erg/s/gr] generated in gamma and positorons
%               from both Ni and Co decay.
%          'Qgamma' - Total energy [erg/s/gr] generated in gammas.
%          'Qpos'   - Total energy [erg/s/gr] generated in positrons.
%          'Q_Ni'   - Total energy from Ni decay.
%          'Q_Co'   - Total energy from Co decay.
%          'f_dep'  - Fraction of energy deposited from the gammas.
%          'Edep'   - Total radiated energy (Qgamma.*f_dep + Qpos).
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Aug 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference : Kulkarni (2006), Eq. 14, 44.
%             Sharon & Kushnir (2020)
% Example: [Res]=AstroUtil.supernova.nickel56_decay([1;2]);
% Reliable: 2
%--------------------------------------------------------------------------

if nargin<3
    n = 3;
    if nargin<2
        t0 = 30;
    end
end


Msun   = 1.9891e33;
Tau_Ni = 8.76;
Tau_Co = 111.4;

L_Ni = 1./Tau_Ni;
L_Co = 1./Tau_Co;
E_Ni = 3.9e10.*exp(-L_Ni.*Time);
E_Co = 7e9.*(exp(-L_Co.*Time) - exp(-L_Ni.*Time));

%E = E_Ni + E_Co;

% new code



Qgamma = (6.54.*exp(-Time./Tau_Ni) + 1.38.*exp(-Time./Tau_Co)).*1e43./Msun;
Qpos   = 4.64e41.*(exp(-Time./Tau_Co) - exp(-Time./Tau_Ni))./Msun;
E      = Qgamma + Qpos;

f_dep = 1./((1 + (Time./t0).^n).^(2./n));


Res.E      = E;
Res.Qgamma = Qgamma;
Res.Qpos   = Qpos;
Res.Q_Ni   = E_Ni;
Res.Q_Co   = E_Co;
Res.f_dep  = f_dep;
Res.Edep   = Qgamma.*f_dep + Qpos;
