function [BestPar,Res]=fit_broken_powlaw(Data,varargin)
% SHORT DESCRIPTION HERE
% Package: Util.fit
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [BestPar,BestChi2]=Util.fit.fit_broken_powlaw(LL(:,1),LL(:,2),LL(:,3));
% Reliable: 
%--------------------------------------------------------------------------


DefV.tbreak               = [0.5 6];
DefV.IndPL                = [0 -1.3 -3.3];
DefV.L0                   = 1e42;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);



Options = optimset('MaxFunEvals',1e6);

Par1 = InPar.tbreak;
Par2 = InPar.IndPL;
Par3 = InPar.L0;

t = Data(:,1);
% L=broken_powlaw(t,tbreak,IndPL,L0)
switch numel(InPar.tbreak)
    case 1
        FunChi2 = @(Par) sum((Util.fun.broken_powlaw(t,Par(1),Par(2:3),Par(4)) - Data(:,2)).^2./(Data(:,3).^2));
    case 2
        FunChi2 = @(Par) sum((Util.fun.broken_powlaw(t,Par(1:2),Par(3:5),Par(6)) - Data(:,2)).^2./(Data(:,3).^2));
    case 3
        FunChi2 = @(Par) sum((Util.fun.broken_powlaw(t,Par(1:3),Par(4:7),Par(8)) - Data(:,2)).^2./(Data(:,3).^2));
    otherwise
        error('1-3 breaks are supported');
end
%Par0 = [6 -1 -3 1e42];
Par0 = [Par1(:);Par2(:);Par3(:)]';


[BestPar,BestChi2,Exit,Out] = fminsearch(FunChi2,Par0,Options);

Res.BestChi2 = BestChi2;
Res.Npar     = numel(BestPar);
Res.Ndof     = numel(t) - Res.Npar;

Res.Hessian    = Util.fit.calc_hessian(FunChi2,BestPar);
Res.BestParErr = sqrt(diag(inv(0.5.*Res.Hessian)));
