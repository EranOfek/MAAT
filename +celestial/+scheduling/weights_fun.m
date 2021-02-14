function W=weights_fun(Fun,t,Pars)
%
%   'exp' - A + B*exp(-t/C)
%   'fermi' - A + B/(1+ exp(-(t-C)/D))

% Example: W=celestial.scheduling.weights_fun('fermiexp',(0:0.01:10).',[1,0,1.5,1,0.1,1,0.5,1])


switch lower(Fun)
    case 'exp'
        W = Pars(1) + Pars(2).*exp(-t./Pars(3));
    case 'fermi'
        W = Pars(1) + Pars(2)./(1 + exp(-(t-Pars(2))./Pars(4)));
    case 'fermiexp'
        W = zeros(size(t));
        W(t<Pars(1))  = Pars(2) + Pars(3)./(1 + exp(-(t(t<Pars(1))-Pars(4))./Pars(5)));
        W(t>=Pars(1)) = Pars(6) + Pars(7).*exp(-t(t>=Pars(1))./Pars(8));
    otherwise
        error('Unknown Fun option');
end