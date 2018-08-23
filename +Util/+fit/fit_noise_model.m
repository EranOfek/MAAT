function Res=fit_noise_model(Resid,Mag)
% SHORT DESCRIPTION HERE
% Package: Util
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

ZP = 27;
StepMag = 0.1;

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Flux = 10.^(-0.4.*(Mag-ZP));

MagVec = (min(Mag):StepMag:max(Mag))';
Nvec   = numel(MagVec);

RR  = nan(Nvec,1);
Par = nan(2,1);
MinRes = Inf;
N = numel(Resid);
H = [ones(N,1), 1./sqrt(Flux)];

for Ivec=1:1:Nvec
    
    Flag = Mag>MagVec(Ivec);
    
    N = numel(Resid(Flag));
    if (N>3)
        Par = H(Flag,:)\Resid(Flag);
    
        Residuals = H(Flag,:)*Par - Resid(Flag);
        RR(Ivec)  = std(Residuals);
        
        if (RR(Ivec)<MinRes)
            BestPar = Par;
        end
    end
    
end

[~,I] = min(RR);

Res.Par    = BestPar;
Res.MagCut = MagVec(I);
