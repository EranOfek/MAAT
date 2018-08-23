function [ResPMP]=fit_pm_parallax(JD,OffRA,OffDec,varargin)
%------------------------------------------------------------------------------
% fit_pm_parallax function                                              FitFun
% Description: Fit low-accuracy parallax and proper motion to a set of
%              of celestial positions.
% Input  : - Vector of JDs on which the position measurments were obtained.
%          - One or two columns matrix in which the first column is the RA
%            offset from a reference position (any consistent units),
%            and the optional second column is for the errors.
%            Note that by the units must be in arcs (e.g., radians,
%            arcseconds), and NOT in arcs of time (e.g., second of times).
%          - One or two columns matrix in which the first column is the Dec
%            offset from a reference position (any consistent units),
%            and the optional second column is for the errors.
%          * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'RA'  - Mean RA or Vector of J2000.0 R.A.
%                    (see convertdms.m for options).
%                    This is used for calculating the parallax direction.
%                    If not provided or empty than assumes that the the
%                    second input argument is given in RA units.
%                    Default is empty.
%            'Dec' - Mean Dec or Vector of J2000.0 Dec.
%                    (see convertdms.m for options).
%                    This is used for calculating the parallax direction.
%                    If not provided or empty than assumes that the the
%                    second input argument is given in Dec units.
%                    Default is empty.
%            'CE'  - Vector of addiative cosmic errors to try in order
%                    to set the \chi^2/dof to 1.
%                    Default is [0:0.01:0.05].'.
%            'ErrRA' - Set error in RA offset. Default is 1.
%            'ErrDec'- Set error in Dec offset. Default is 1.
%            'Plx'   - Fit parallax {true|false}. Default is true.
%            'ModelJD' - If a list of JD is provided than will add to
%                    ouput the model calculated at these JD.
%                    Default is empty matrix.
% Output : - Best fit parameters:
%            [best OffsetRA (J2000),
%             best PM RA [in time units]
%             best parralax (from RA)
%             best OffsetDec (J2000),
%             best PM Dec 
%             best parralax (from Dec)]
%          - Errors in best fit parameters.
%          - \chi^2
%          - d.o.f.
%          - Vector of residuals [RA followed by Dec].
% Tested : Matlab 7.6
%     By : Eran O. Ofek                 September 2008
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------
J2000= 2451545.0;   % (J2000.0)
JYear= 365.25;      % Julian year

DefV.RA     = [];
DefV.Dec    = [];
DefV.CE     = [0:0.01:1].';
DefV.ErrRA  = 1;
DefV.ErrDec = 1;
DefV.Plx    = true;
DefV.ModelJD= [];
InPar = set_varargin_keyval(DefV,'y','use',varargin{:}); 

if (isempty(InPar.RA)),
    RA  = OffRA;
end
if (isempty(InPar.Dec)),
    Dec = OffDec;
end


FlagNN = ~isnan(JD) & ~isnan(OffRA(:,1)) & ~isnan(OffDec(:,1)); % & ~isnan(RA) & ~isnan(Dec);
JD     = JD(FlagNN);
OffRA  = OffRA(FlagNN);
OffDec = OffDec(FlagNN);
%RA     = RA(FlagNN);
%Dec    = Dec(FlagNN);

Nobs = numel(JD);

if (size(OffRA,2)==1),
    OffRA  = [OffRA, ones(Nobs,1).*InPar.ErrRA];
end
if (size(OffDec,2)==1),
    OffDec = [OffDec, ones(Nobs,1).*InPar.ErrDec];
end



RA   = convertdms(InPar.RA,'gH','r');
Dec  = convertdms(InPar.Dec,'gD','R');

% J2000.0 Equatorial position of the Earth [au]
if (InPar.Plx),
   [Coo,Vel] = calc_vsop87(JD, 'Earth', 'e', 'E');
   X = Coo(1,:).';
   Y = Coo(2,:).';
   Z = Coo(3,:).';
end

% alpha1 ~ alpha + pm_alpha.*(T-T0)
%                + pi.*(X.*sin(alpha1) - Y.*cos(./(15.*cos(delta))


% Mean position
RA_t  = mean(InPar.RA);
Dec_t = mean(InPar.Dec);

% design matrix for RA + PM in RA
Ha = [ones(Nobs,1), (JD-J2000)./JYear];
Ya = OffRA(:,1);
Yae= OffRA(:,2);

% design matrix for Dec
Hd = [ones(Nobs,1), (JD-J2000)./JYear];
Yd = OffDec(:,1);
Yde= OffDec(:,2);

if (InPar.Plx),
    % add Parallax to design matrix
    %Hp = [(X.*sin(RA_t) - Y.*cos(RA_t))./(15.*cos(Dec_t)); (X.*cos(RA_t).*sin(Dec_t) + Y.*sin(RA_t).*sin(Dec_t) - Z.*cos(Dec_t))];
    Hp = [(X.*sin(RA_t) - Y.*cos(RA_t)); (X.*cos(RA_t).*sin(Dec_t) + Y.*sin(RA_t).*sin(Dec_t) - Z.*cos(Dec_t))];
else
    Hp = zeros(2.*Nobs,0);
end

% design matrix
H = [[Ha, zeros(Nobs,2); zeros(Nobs,2), Hd], Hp];
Y = [Ya; Yd];
Ye= [Yae;Yde];
Na = length(Ya);
Nd = length(Yd);

Nce  = numel(InPar.CE);
Chi2 = zeros(Nce,1);
for Ice=1:1:Nce,
    Yvar       = Ye.^2 + InPar.CE(Ice).^2;
    [Par]      = lscov(H,Y,1./Yvar);
    Resid      = Y - H*Par;
    Chi2(Ice)  = sum(Resid.^2./Yvar);
end
Dof          = Nobs.*2 - 5;
ListZ = find_local_zeros(InPar.CE,Chi2./Dof - 1);
if (isempty(ListZ)),
    % set cosmic error to 0
    if (Chi2./Dof>2),
        CE = max(InPar.CE);
    else
       CE = 0;
    end
else
    CE = min(ListZ(:,1));
end
Yvar          = Ye.^2 + CE.^2;
[Par,ParErr]  = lscov(H,Y,1./Yvar);

Resid         = Y - H*Par;
Chi2          = sum(Resid.^2./Yvar);

%[~,~,Flag]=clip_resid(Resid,varargin);



ResPMP.Par    = Par;
ResPMP.ParErr = ParErr;
ResPMP.Cov    = inv(H'*inv(diag(Yvar))*H);
ResPMP.CE     = CE;
ResPMP.Chi2   = Chi2;
ResPMP.Dof    = Dof;
ResPMP.Nobs   = Nobs;
ResPMP.Resid  = Resid;
ResPMP.ResidX = Resid(1:Na);
ResPMP.ResidY = Resid(Na+1:end);

ResPMP.ModelJD = InPar.ModelJD;


% Return design matrix for plotting
if (isempty(InPar.ModelJD)),
    ResPMP.ModelH  = [];
    ResPMP.ModelX  = [];
    ResPMP.ModelY  = [];
else
    % J2000.0 Equatorial position of the Earth [au]
    if (InPar.Plx),
        [Coo,Vel] = calc_vsop87(InPar.ModelJD, 'Earth', 'e', 'E');
        X = Coo(1,:).';
        Y = Coo(2,:).';
        Z = Coo(3,:).';
    end
    Nmodel = numel(InPar.ModelJD);
    % design matrix for RA + PM in RA
    Ha = [ones(Nmodel,1), (InPar.ModelJD-J2000)./JYear];

    % design matrix for Dec
    Hd = [ones(Nmodel,1), (InPar.ModelJD-J2000)./JYear];

    if (InPar.Plx),
        % add Parallax to design matrix
        %Hp = [(X.*sin(RA_t) - Y.*cos(RA_t))./(15.*cos(Dec_t)); (X.*cos(RA_t).*sin(Dec_t) + Y.*sin(RA_t).*sin(Dec_t) - Z.*cos(Dec_t))];
        Hp = [(X.*sin(RA_t) - Y.*cos(RA_t)); (X.*cos(RA_t).*sin(Dec_t) + Y.*sin(RA_t).*sin(Dec_t) - Z.*cos(Dec_t))];
    else
        Hp = zeros(2.*Nmodel,0);
    end

    % design matrix
    ResPMP.ModelH = [[Ha, zeros(Nmodel,2); zeros(Nmodel,2), Hd], Hp];
    ModelY        = ResPMP.ModelH*ResPMP.Par;
    ResPMP.ModelX = ModelY(1:Nmodel);
    ResPMP.ModelY = ModelY(Nmodel+1:end);
end
