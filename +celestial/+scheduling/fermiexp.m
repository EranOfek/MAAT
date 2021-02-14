function W=fermiexp(t,t0,Decay,Soft,BaseW,ExtraW)
% Fermi rise - Exp decay weight function
% Package: +celestial.scheduling
% Input  : - time
%          * t0, DecayExp, SoftFermi, BaseW, ExtraW
% Output : - Weights
% Example: W=celestial.scheduling.fermiexp(t,1,1,0.05,1,0.5);

W = zeros(size(t));
W(t<t0)  = (BaseW + ExtraW)./(1 + exp(-(t(t<t0)-t0)./Soft));
W(t>=t0) = BaseW + ExtraW.*exp(-(t(t>=t0)-t0)./Decay);
        
        
        