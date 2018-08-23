function L=broken_powlaw(t,tbreak,IndPL,L0)
% Generate a broken power-law time series.
% Package: Util.fun
% Description: Generate a broken power-law time series with arbitrary
%              number of breaks.
% Input  : - Vector of times.
%          - Vector of break times.
%          - Vector of power-law indices (larger by 1 than the number of
%            breaks).
%          - Normalization at t=1.
% Output : - Vector of time series.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: t=(0.1:0.1:20)'; L=Util.fun.broken_powlaw(t,[1 6],[-0.1 -1.1 -3],1e42)
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==2)
    All = tbreak;
    L0  = All(end);
    N   = numel(All);
    Nslope = N./2;
    tbreak = All(1:Nslope-1);
    IndPL  = All(Nslope:end-1);
end
    

Nbreak = numel(tbreak);

if (Nbreak+1)~=numel(IndPL)
    error('Number of tbreak should be smaller by 1 than number of IndPL');
end

L = L0.*t.^IndPL(1);
Li = L0;
for Ibreak=1:1:Nbreak
    It = find(t>tbreak(Ibreak));
    ti = tbreak(Ibreak);
    Li = Li.*ti.^(IndPL(Ibreak))./(ti.^IndPL(Ibreak+1));
    L(It) = Li.*t(It).^IndPL(Ibreak+1);
end

