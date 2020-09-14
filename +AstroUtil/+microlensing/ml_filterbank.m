function [F,Pars]=ml_filterbank(varargin)
% Generate microlensing template bank
% Package: AstroUtil
% Description: Generate microlensing template bank
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            't' - Vector of times. Default is (-100:1:100).'.
%            'VecBeta' - Vector of minmal impact parameters.
%                        Default is (0.1:0.05:1).^2.'.
%            'VecV'    - Vector of velocty in Einstein radius per day.
%                        Default is logspace(-3,-0.5,25).'.
%            'OutMag'  - Output is in magnitude scaling. Default is true.
% Output : - Matrix of templates bank.
%            Each column correspond to one template.
%          - Structure containing value of parameters in each column of the
%            template bank.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [F,ParsML]=AstroUtil.microlensing.ml_filterbank;
% Reliable: 
%--------------------------------------------------------------------------



DefV.t                   = (-100:1:100).';
DefV.VecBeta             = (0.1:0.05:1).^2.';
DefV.VecV                = logspace(-3,-0.5,25).';
DefV.OutMag              = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


NV = numel(InPar.VecV);
NB = numel(InPar.VecBeta);
Nt = numel(InPar.t);

F = zeros(Nt,NV.*NB);
Pars.VecBeta = zeros(1,NV.*NB);
Pars.VecV    = zeros(1,NV.*NB);

K = 0;
for Iv=1:1:NV
    for Ib=1:1:NB
        K = K + 1;
        [Mag,Res]=AstroUtil.microlensing.microlens_ps([0 InPar.VecBeta(Ib) InPar.VecV(Iv) 1 0],InPar.t);
        if (InPar.OutMag)
            F(:,K) = 2.5.*log10(Res.Mu);
        else
            F(:,K) = Res.Mu - 1;
        end
        Pars.VecBeta(K) = InPar.VecBeta(Ib);
        Pars.VecV(K)    = InPar.VecV(Iv);
    end
end