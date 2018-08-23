function [Res,Fun]=ps_lens(varargin)
% Calculate deflection, magnification and time delay for point mass lens
% Package: AstroUtil
% Description: Calculate deflection, magnification and time delay for point
%              mass lens.
% Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Mass' - Lens mass. Default is 1.
%            'MassUnits' - Lens mass units. Default is 'SunM'
%            'Dl'   - Lens distance or redshift. Default is 5000;
%            'Ds'   - Source distance or redshift. Default is 1e4;
%            'DistUnits' - Distance units. Default is 'pc'.
%            'Beta' - Impact parameter. Default is logspace(-3,3,7)'.
%            'BetaUnits' - Impact parameter units. Default is 'ThetaE'.
%            'BetaMin' - Minimum impact parameter. Default is 0.
%                       If provided than the lens source distance is
%                       calculated from sqrt(Beta^2 + BetaMin^2). BetaMin
%                       has the units specified in BetaUnits.
%            'OutUnits' - Units of angular output. Default is 'arcsec'.
% Output : - A structure containing the microlensing properties.
%          - A structure of useful functions.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=AstroUtil.microlensing.ps_lens
% Reliable: 2
%--------------------------------------------------------------------------


DefV.Mass                 = 1;
DefV.MassUnits            = 'SunM';
DefV.Dl                   = 5000;
DefV.Ds                   = 1e4;
DefV.DistUnits            = 'pc'; % 'z' | 'pc' | 'cm'
DefV.Beta                 = logspace(-3,3,7)';
DefV.BetaUnits            = 'ThetaE'; % 'rad' | 'arcsec' | 'mas' | 'ThetaE'  
DefV.BetaMin              = 0;
DefV.OutUnits             = 'arcsec'; % 'rad' | 'arcsec' | ...
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% unit conversion
G  = constant.G;
c  = constant.c;
pc = constant.pc; 
switch lower(InPar.DistUnits)
    case 'z'
        Dl  = AstroUtil.cosmo.ad_dist(InPar.Dl);
        Ds  = AstroUtil.cosmo.ad_dist(InPar.Ds);
        Dls = AstroUtil.cosmo.ad_dist([InPar.Dl,InPar.Ds]);
        
        % output is pc
    otherwise
        Dl  = InPar.Dl;
        Ds  = InPar.Ds;
        Dls = InPar.Ds - InPar.Dl;
        switch lower(InPar.DistUnits)
            case 'pc'
                % already in pc
            case 'cm'
                % convert cm to pc
                Dl  = Dl./pc;
                Ds  = Ds./pc;
                Dls = Dls./pc;
            case 'kpc'
                % convert kpc to pc
                Dl  = Dl.*1000;
                Ds  = Ds.*1000;
                Dls = Dls.*1000;
            case 'mpc'
                % convert Mpc to pc
                Dl  = Dl.*1e6;
                Ds  = Ds.*1e6;
                Dls = Dls.*1e6;    
            otherwise
                error('Unknwon InUnits option');
        end
end

M  = convert.mass(InPar.MassUnits,'gr',InPar.Mass);

FunThetaE = @(M,Dl,Ds,Dls) sqrt( (4.*G.*M./(c.^2)) .* Dls./(Dl.*Ds.*pc) );
TE = FunThetaE(M,Dl,Ds,Dls);  % Einstein radius in radians

switch lower(InPar.BetaUnits)
    case 'thetae'
        % convert beta from ThetaE to radians
        Beta    = InPar.Beta.*TE;
        BetaMin = InPar.BetaMin.*TE;
    otherwise
        % convert beta from angular units to radians
        Beta    = convert.angular(InPar.BetaUnits,'rad',InPar.Beta);
        BetaMin = convert.angular(InPar.BetaUnits,'rad',InPar.BetaMin);
end

% add BetaMin
Beta = sqrt(Beta.^2 + BetaMin.^2);  % radians

FunThetaP = @(beta,ThetaE) 0.5.*(beta + sqrt(beta.^2 +4.*ThetaE.^2));
FunThetaM = @(beta,ThetaE) 0.5.*(beta - sqrt(beta.^2 +4.*ThetaE.^2)); 
% note magnification may be negative (abs value should be taken)
FunMuP = @(beta,ThetaE) 1./(1 - (ThetaE./FunThetaP(beta,ThetaE)).^4 );
FunMuM = @(beta,ThetaE) 1./(1 - (ThetaE./FunThetaM(beta,ThetaE)).^4 );

ThetaP = FunThetaP(Beta,TE);
ThetaM = FunThetaM(Beta,TE);
MuP    = FunMuP(Beta,TE);
MuM    = FunMuM(Beta,TE);

ShiftReLens   = (abs(MuP).*ThetaP + abs(MuM).*ThetaM)./(abs(MuP) + abs(MuM));
ShiftReSource = ShiftReLens - Beta;

% convert angular values to output units
ConvAng = convert.angular('rad',InPar.OutUnits);

Res.ER      = TE.*ConvAng;
Res.ThetaP  = ThetaP.*ConvAng;
Res.ThetaM  = ThetaM.*ConvAng;
Res.MuP     = MuP;
Res.MuM     = MuM;
Res.MuTot   = abs(MuP) + abs(MuM);
Res.ShiftReLens   = ShiftReLens.*ConvAng;
Res.ShiftReSource = ShiftReSource.*ConvAng;
Res.Beta          = Beta.*ConvAng;

% potential
Phi1    = TE.^2.*log(abs(ThetaP));
Phi2    = TE.^2.*log(abs(ThetaM));

% time delay
D   = Dls./(Dl.*Ds.*pc);
switch lower(InPar.DistUnits)
    case 'z'
        % note: InPar.Dl contains redshift
        Res.Delay1 = (1 + InPar.Dl).*(0.5.*(ThetaP - Beta).^2 - Phi1)./(c.*D);
        Res.Delay2 = (1 + InPar.Dl).*(0.5.*(ThetaM - Beta).^2 - Phi2)./(c.*D);
    otherwise
        % assume not cosmological time dilation
        %Res.Delay1 = (0.5.*(ThetaP - Beta).^2 - Phi1)./(c.*D);
        %Res.Delay2 = (0.5.*(ThetaM - Beta).^2 - Phi2)./(c.*D);
        A = -2.*G.*M./(c.^3);
        Res.Delay1 = A .* log(1 - cos(Beta));
        Res.Delay2 = A .* log(1 - cos(Beta));
        
%         Res.dDelay1_db1 = A.* (-sin(Beta))./(cos(Beta)-1);
%         Res.dDelay1_db2 = A./(cos(Beta)-1);
%         Res.dDelay1_db3 = A.*sin(Beta)./((cos(Beta)-1).^2);
%         [dbdt1,dbdt2,dbdt3] = Util.fun.numerical_diff(Res.Beta);
%         Res.dDelay1_dt1 = Res.dDelay1_db1.*dbdt1;
%         Res.dDelay1_dt2 = Res.dDelay1_db2.*dbdt2;
%         Res.dDelay1_dt3 = Res.dDelay1_db3.*dbdt3;
end

if (nargout>1)
    Fun.FunThetaP = FunThetaP;
    Fun.FunThetaM = FunThetaM;
    Fun.FunMuP    = FunMuP;
    Fun.FunMuM    = FunMuM;
end
