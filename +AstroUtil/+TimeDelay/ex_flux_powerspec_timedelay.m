function [EP,Omega]=ex_flux_powerspec_timedelay(varargin)
% 
% Package: AstroUtil
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [EP,Omega]=AstroUtil.TimeDelay.ex_flux_powerspec_timedelay;
% Reliable: 
%--------------------------------------------------------------------------


DefV.Freq                 = (0:0.001:1).';
DefV.Tau                  = 20;
DefV.A1                   = 1;
DefV.A2                   = 0.5;
DefV.gamma                = 3;
DefV.Err                  = 0.01;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Omega = 2.*pi.*InPar.Freq;
EP = (InPar.A1.^2 + InPar.A2.^2 + 2.*InPar.A1.*InPar.A2.*cos(Omega.*InPar.Tau))./InPar.Freq.^InPar.gamma + InPar.Err.^2;
