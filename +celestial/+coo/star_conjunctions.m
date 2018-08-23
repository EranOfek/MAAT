function Res=star_conjunctions(Target1,Target2,varargin)
% Calculate conjuctions between stars given their proper motion.
% Package: celestial.coo
% Description: Given a star with its coordinates, proper motion and
%              optionally parallax and radial velocity, and a list of
%              multiple stars at the same sky region, calculate possible
%              conjunctions between the star in the first list and all the
%              stars in the second list. For each possible conjunction,
%              calculate also the microlensing properties of the event
%              including the astrometric microlensing.
% Input  : - A matrix containing one line, or a an AstCat object containing
%            a single star.
%          - A matrix in which each line correspond to a star, or an AstCat
%            object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            ''
% Output : - A structure containing all possible conjunctions.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: celestial.coo.proper_motion_nearest_approach
% Reliable: 2
%--------------------------------------------------------------------------

RAD        = 180./pi;
ARCSEC_DEG = 3600;
MAS_AS     = 1000;
MAS_DEG    = ARCSEC_DEG.*MAS_AS;
JYear      = 365.25;

CatField     = AstCat.CatField;
ColCellField = AstCat.ColCellField;


DefV.MinDistThreshold     = 30;      % If min dist below this then calc microlensing properties
DefV.Mass                 = 1;       % [solar mass]
DefV.MinAssymDef          = 100e-6;  % [arcsec]
DefV.MinDistEpoch0        = 1.0;       % [arcsec]

%
DefV.ColCell1             = {};
DefV.ColCell2             = {};
%
DefV.CooUnits             = 'rad';  % errors must be in mas
DefV.VecRelYear           = (-10:0.005:10)';
DefV.VecAbsYear           = [];
DefV.MinDistYearRange     = [2000 2070];
DefV.Dist_NoPlx           = 10;   % kpc
DefV.PlxNsigma            = 2;
DefV.DefPlxRelErr         = 0.3;
DefV.DistUnits            = 'kpc';
DefV.LensMass             = 1.4;   % solar mass
%
DefV.ColRA_1              = 'RA';
DefV.ColDec_1             = 'Dec';
DefV.ColPMRA_1            = 'PMRA';
DefV.ColPMDec_1           = 'PMDec';
DefV.ColEpoch_1           = 'PosEpoch';
DefV.ColErrRA_1           = [];
DefV.ColErrDec_1          = [];
DefV.ColErrPMRA_1         = [];
DefV.ColErrPMDec_1        = [];
DefV.ColPlx_1             = [];
DefV.ColErrPlx_1          = [];
DefV.ColDist_1            = 'Dist';  % kpc
%
DefV.ColRA_2              = 'RA';
DefV.ColDec_2             = 'Dec';
DefV.ColPMRA_2            = 'PMRA';
DefV.ColPMDec_2           = 'PMDec';
DefV.ColEpoch_2           = 'Epoch';
DefV.ColErrRA_2           = 'ErrRA';
DefV.ColErrDec_2          = 'ErrDec';
DefV.ColErrPMRA_2         = 'ErrPMRA';
DefV.ColErrPMDec_2        = 'ErrPMDec';
DefV.ColPlx_2             = 'Plx';
DefV.ColErrPlx_2          = 'ErrPlx';
DefV.ColDist_2            = [];
%
DefV.UnitsEpoch_1         = 'MJD';
DefV.UnitsEpoch_2         = 'J';


InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% convert Target1/Target2 to AstCat objects
if (~AstCat.isastcat(Target1))
    Cat(1) = AstCat;
    Cat(1).(CatField) = Target1;
    if (isempty(InPar.ColCell1))
        error('Cat1 is an AstCat object - must provide ColCell1');
    end
    Cat(1).(ColCellField) = InPar.ColCell1;
    Cat(1) = colcell2col(Cat(1));
else
    Cat(1) = Target1;
end

if (~AstCat.isastcat(Target2))
    Cat(2) = AstCat;
    Cat(2).(CatField) = Target2;
    if (isempty(InPar.ColCell2))
        error('Cat2 is an AstCat object - must provide ColCell2');
    end
    Cat(2).(ColCellField) = InPar.ColCell2;
    Cat(2) = colcell2col(Cat(2));
else
    Cat(2) = Target2;
end

if (size(Cat(1).(CatField),1)>1)
    error('Cat1 should contain a single source');
end

% extract relevant field from Cat1/Cat2
Fields = {'RA','Dec','PMRA','PMDec','ErrRA','ErrDec','ErrPMRA','ErrPMDec','Epoch','Plx','ErrPlx','Dist'};
Nf     = numel(Fields);

for Icat=1:1:2
    Ns = size(Cat(Icat).(CatField),1);
    for If=1:1:Nf
        %[Icat,If]
        [FieldName,ColName] = get_colname(Fields{If},Icat);
        
        if (isempty(InPar.(ColName)))
            T(Icat).(FieldName) = zeros(Ns,1);
        else
            T(Icat).(FieldName) = col_get(Cat(Icat),InPar.(ColName));
        end
    end
end

% convert epoch to Julian Years
T(1).Epoch = convert.time(T(1).Epoch,InPar.UnitsEpoch_1,'J');
T(2).Epoch = convert.time(T(2).Epoch,InPar.UnitsEpoch_2,'J');

T(1).RA  = convert.angular(InPar.CooUnits,'rad',T(1).RA);
T(2).RA  = convert.angular(InPar.CooUnits,'rad',T(2).RA);
T(1).Dec = convert.angular(InPar.CooUnits,'rad',T(1).Dec);
T(2).Dec = convert.angular(InPar.CooUnits,'rad',T(2).Dec);

% treat distance and parallax
for It=1:1:2
    % convert ditsnace to kpc
    switch lower(InPar.DistUnits)
        case 'kpc'
            % do nothing
        case 'pc'
            % pc to kpc
            T(It).Dist = T(It).Dist./1000;
        otherwise
            error('Unknown DistUnits option');
    end

    Ns = numel(T(It).RA);
    for Is=1:1:Ns
        if (T(It).Plx(Is)<=(InPar.PlxNsigma.*T(It).ErrPlx(Is)))
            if (T(It).Dist>0)
                % T(1).Dist is defined, but no parallax
                T(It).Plx(Is)    = 1./(T(It).Dist(Is));  % mas
                T(It).ErrPlx(Is) = T(It).Plx(Is).*InPar.DefPlxRelErr;
            else
                % no parallax and no distance
                T(It).Dist(Is)   = InPar.Dist_NoPlx;  % kpc
                T(It).Plx(Is)    = 1./T(It).Dist(Is); % mas
                T(It).ErrPlx(Is) = T(It).Plx(Is).*InPar.DefPlxRelErr;
            end
        else
            % parallax is defined
            if (T(It).Dist(Is)<=0)
                % no distance
                T(It).Dist(Is) = 1./T(It).Plx(Is);  % kpc
            end
        end
    end
end
    
    

% convert to mas/yr units:
RA2mas  = RAD.*MAS_DEG.*cos(T(1).Dec);
Dec2mas = RAD.*MAS_DEG;

DistEpoch0 = celestial.coo.sphere_dist_fast(T(1).RA,T(1).Dec,T(2).RA,T(2).Dec).*RAD.*ARCSEC_DEG;   % [arcsec]

[MinT,MinDist]=Util.Geom.traj_mindist(T(1).RA.*RA2mas, T(1).PMRA, T(1).Dec.*Dec2mas, T(1).PMDec, T(1).Epoch,...
                                      T(2).RA.*RA2mas, T(2).PMRA, T(2).Dec.*Dec2mas, T(2).PMDec, T(2).Epoch);


MinDist        = MinDist./MAS_AS;  % arcsec
MinDist_rad    = MinDist./(RAD.*ARCSEC_DEG);
Res.AllMinT    = MinT;
Res.AllMinDist = MinDist;
% if (isempty(Res.AllMinDist))
%     Res.MinMinDist = NaN;
% else
%     Res.MinMinDist = min(MinDist);
% end

% Estimate asymptotic light deflection:
% assuming beta>>ThetaE and Ds=infinity
Dl = 1000./T(1).Plx;  % [pc]
Assym_def = 4.*constant.G.*constant.SunM.*InPar.Mass./(constant.c.^2 .* Dl.*constant.pc .* MinDist_rad);  % [rad]
Assym_def = Assym_def.*RAD.*ARCSEC_DEG;  % [arcsec]

Ifound = find(DistEpoch0>InPar.MinDistEpoch0 & Assym_def>InPar.MinAssymDef & MinDist<(InPar.MinDistThreshold) & MinT>InPar.MinDistYearRange(1) & MinT<InPar.MinDistYearRange(2));

%Res = Util.struct.struct_def({'Found','Nfound','Ifound','MinT','MinDist','Cat1','Cat2','Sol'},1,1);

if (isempty(Ifound))
    Res.Found    = false;
    Res.Nfound   = 0;
    Res.Ifound   = [];
    Res.Cat1     = Cat(1);
    Res.Cat2     = Cat(2);
    Res.ConjTime = MinT;
    Res.ConjDist = MinDist;
    Res.Cand     = [];
    Res.bestDist       = NaN;
    Res.bestDrange     = NaN;
    Res.bestDeflection = NaN;
    Res.bestMagnification = NaN;
    
%     Res.MinMinDeflection = NaN;
%     Res.MaxShiftDeflection = NaN;
%     Res.rangeDmax  = NaN;
%     Res.MinMinDist = NaN;
else
    Res.Found        = true;
    Res.Nfound       = numel(Ifound);
    Res.Ifound       = Ifound;
    Res.Cat1         = Cat(1);
    Res.Cat2         = Cat(2);
    Res.ConjTime     = MinT(Ifound);
    Res.ConjDist     = MinDist(Ifound);
    
    for Iff=1:1:Res.Nfound
        % for each candidate calc detailes
        If2 = Ifound(Iff);
        
        if (isempty(InPar.VecAbsYear))
            EpochOut = Res.ConjTime(Iff) + InPar.VecRelYear;
        else
            EpochOut = InPar.VecAbsYear;
        end
        
        % Note epoch units should be days
        [RA1,Dec1] = celestial.coo.proper_motion(EpochOut.*JYear,T(1).Epoch.*JYear,      T(1).Epoch.*JYear,      T(1).RA,      T(1).Dec,      T(1).PMRA,      T(1).PMDec,      T(1).Plx);
        [RA2,Dec2] = celestial.coo.proper_motion(EpochOut.*JYear,T(2).Epoch(If2).*JYear, T(2).Epoch(If2).*JYear, T(2).RA(If2), T(2).Dec(If2), T(2).PMRA(If2), T(2).PMDec(If2), T(2).Plx(If2));
        
        Res.Cand(Iff).RA1  = RA1;
        Res.Cand(Iff).Dec1 = Dec1;
        Res.Cand(Iff).RA2  = RA2;
        Res.Cand(Iff).Dec2 = Dec2;
        
        [D,PA] = celestial.coo.sphere_dist_fast(RA1,Dec1,RA2,Dec2);

        Res.Cand(Iff).EpochOut   = EpochOut;        % JYear
        Res.Cand(Iff).DistT      = D.*RAD.*ARCSEC_DEG;  % arcsec
        Res.Cand(Iff).PA         = PA.*RAD;    % deg
        
        Res.Cand(Iff).rangeDist  = range(D).*RAD.*ARCSEC_DEG;  % arcsec
        [Res.Cand(Iff).minDist,Imin]  = min(Res.Cand(Iff).DistT); % arcsec
       
        
        Res.Cand(Iff).ErrRA    = sqrt( T(1).ErrRA.^2 + ((T(1).Epoch - EpochOut).*T(1).ErrPMRA).^2 + ...
                                           T(2).ErrRA(If2).^2 + ((T(2).Epoch(If2) - EpochOut).*T(2).ErrPMRA(If2)).^2  + 500.^2);   % mas
        Res.Cand(Iff).ErrDec   = sqrt( T(1).ErrDec.^2 + ((T(1).Epoch - EpochOut).*T(1).ErrPMDec).^2 + ...
                                           T(2).ErrDec(If2).^2 + ((T(2).Epoch(If2) - EpochOut).*T(2).ErrPMDec(If2)).^2 + 500.^2);        % mas                       
            
        
        Res.Cand(Iff).ErrD      = sqrt(Res.Cand(Iff).ErrRA.^2 + Res.Cand(Iff).ErrDec.^2);  % mas
        Res.Cand(Iff).ErrD_minD = Res.Cand(Iff).ErrD(Imin);   % mas
        
        Res.Cand(Iff).RelPMRA  = T(1).PMRA  - T(2).PMRA(If2);
        Res.Cand(Iff).RelPMDec = T(1).PMDec - T(2).PMDec(If2);
        Res.Cand(Iff).RelPM    = sqrt(Res.Cand(Iff).RelPMRA.^2 + Res.Cand(Iff).RelPMDec.^2);
        
        Res.Cand(Iff).Err_MinT    = sqrt(Res.Cand(Iff).ErrRA(Imin).^2 + Res.Cand(Iff).ErrDec(Imin).^2);
        Res.Cand(Iff).IndCat2     = If2;
        
        % microlensing signal
        D1 = T(1).Dist.*1000; % pc
        D2 = T(2).Dist(If2).*1000; % pc
        if (isnan(D2))
            D2 = 1e4;
        end
        MinD = min(D1,D2);
        MaxD = max(D1,D2);
        
        D_rad = Res.Cand(Iff).DistT./(RAD.*ARCSEC_DEG);  % angular dist [rad]
        ResM = AstroUtil.microlensing.ps_lens('Mass',InPar.LensMass,'Dl',MinD,'Ds',MaxD,'Beta',D_rad,'BetaUnits','rad');
        
        Res.Cand(Iff).Dist1      = D1;
        Res.Cand(Iff).Dist2      = D2;
        
        Res.Cand(Iff).ER         = ResM.ER;
        Res.Cand(Iff).ML         = ResM;
        
        Res.Cand(Iff).MaxDeflection = max(ResM.ShiftReSource);
        Res.Cand(Iff).MaxMagnification = max(ResM.MuTot);
    end
    
    [Res.bestDeflection,Im] = max([Res.Cand.MaxDeflection]);
    Res.bestDist      = Res.Cand(Im).minDist;
    Res.bestDrange    = Res.Cand(Im).rangeDist;
    Res.bestMagnification = max([Res.Cand.MaxMagnification]);
    
end
     




end
    
function [FieldName,ColName] = get_colname(Prop,Icat)

FieldName = sprintf('%s',Prop);
if (nargout>1)
    ColName = sprintf('Col%s_%d',Prop,Icat);
end

end
