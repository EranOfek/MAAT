function [R,Mz,P]=volumetric_rate_sn(varargin)
% SHORT DESCRIPTION HERE
% Package: telescope.sn
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [d]=telescope.sn.volumetric_rate_sn(varargin)
% Reliable: 
%--------------------------------------------------------------------------



DefV.SFR                  = false;
DefV.Omega                = 1;
DefV.A                    = (2:0.3:8).'; %(0.3:0.3:8)';
DefV.Sigma                = 1;
DefV.SN                   = 5;
DefV.B                    = 1;
DefV.Eta                  = 1;
DefV.W                    = 1;
DefV.t                    = 1;
DefV.L                    = 1e13; %logspace(16,20,10)';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


d = InPar.L.^(1./2) .* InPar.SN.^(-1./2) .* InPar.B.^(-1./4) .* InPar.Sigma.^(-1./2) .* (InPar.Eta.*InPar.A.*InPar.W.*InPar.t).^(1./4);


z = AstroUtil.cosmo.inv_lum_dist(d,'LD');
[~,Vc]=AstroUtil.cosmo.comoving_volume(z);

if (InPar.SFR)
    SFR = AstroUtil.cosmo.sfr(z);

    R = Vc.*InPar.Omega.*SFR;
else
    R = Vc.*InPar.Omega;
end

Mz = mean(z);

 
[min(z),max(z)]
P=Util.fit.fitpow(InPar.A,R,ones(size(R)));

% VecL = logspace(13,22,100)';
% for Il=1:1:numel(VecL)
%     [R,Mz,P]=telescope.sn.volumetric_rate_sn('L',VecL(Il));
%     D(Il,1:2) = [Mz, P(2)];
% end
% plot(D(:,1),D(:,2),'k-','LineWidth',2);
% %hold on; plot(D(:,1),D(:,2)./(1+D(:,1)),'k--','LineWidth',2);
% axis([0 15 0 0.8])
% H =xlabel('z'); H.FontSize=18; H.Interpreter='latex';
% H =ylabel('$\alpha$'); H.FontSize=18; H.Interpreter='latex';
% print alpha_z.eps -depsc2


