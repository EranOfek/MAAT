function [RMS,Factor,Res]=fit_blackbody_spec(TempVec,Spec,FluxUnits,WaveUnits,Dist,Ebv,Ndof)
%----------------------------------------------------------------------------
% fit_blackbody_spec function                                      AstroSpec
% Description: Fit a spectrum to a black-body spectra.
%              OBSOLOTE: Use fit_bb.m instead.
% Input  : - Vector of temperatures to test (e.g., [1000:500:30000]') [K].  
%            If empty matrix the use default. Default is logspace(3,6,30).'
%          - Spectrum to fit [Wave, Flux].
%            Default wavelength is in Angstrem and default flux in
%            erg/cm^/s/A.
%          - Flux units:
%            'F_lambda'  - erg/cm^2/s/A (default).
%            'F_nu'      - erg/cm^2/s/Hz.
%            'mJy'       - mJy.
%            'F_kev'     - erg/cm^2/s/keV
%            'F_ev'      - erg/cm^2/s/eV
%          - Wavelength units:
%            'A'         - Angstrem (default).
%            'nm'        - nanometer.
%            'cm'        - cm.
%            'm'         - meter.
%            'Hz'        - Hz
%            'eV'        - eV
%            'keV'       - keV
%          - Optional distance in pc, default is 10pc.
%          - Extinction correction (E_{B-V}) to apply to the
%            spectrum before the fit assuming the spectrum
%            is extincted at redshift zero. Default is 0.
%          - Number of degrees of freedom to use in the error estimation.
%            Provide the number of uncorrelated points in the spectrum - 2.
%            Default is the number of points in spectrum - 2.
% Output : - Vector of RMS of best fit (as a function of temperature).
%          - Vector of best fit multiplication factor between the spectra
%            divided by a black body spectrum emitted from 1 cm^2.
%            To convert to the source radius [cm]:
%            sqrt(Factor.*Dist.^2), where Dist is the distance to
%            the source in cm.
%          - Structure of best fit parameters, with the following fields:
%            .Ind        - Index of best fit in temperature vector.
%            .BestT      - Interpolated best fit temperature.
%            .BestRMS    - Best fit RMS (not interpolated).
%            .BestFactor - Interpolated best fit factor.
%            .BestRadius - Interpolated best fit radius [cm].
%            If didn't find local minimum then return NaN.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    August 2010     
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%----------------------------------------------------------------------------

Pc = get_constant('pc');

Def.FluxUnits  = 'F_lambda';
Def.WaveUnits  = 'A';
Def.Dist       = 10;   % [pc]
Def.Ebv        = 0;    % [mag]
Def.Ndof       = length(Spec(:,1)) - 2;

if (nargin==2),
   FluxUnits   = Def.FluxUnits;
   WaveUnits   = Def.WaveUnits;
   Dist        = Def.Dist;
   Ebv         = Def.Ebv;
   Ndof        = Def.Ndof;
elseif (nargin==3),
   WaveUnits   = Def.WaveUnits;
   Dist        = Def.Dist;
   Ebv         = Def.Ebv;
   Ndof        = Def.Ndof;
elseif (nargin==4),
   Dist        = Def.Dist;
   Ebv         = Def.Ebv;
   Ndof        = Def.Ndof;
elseif (nargin==5),
   Ebv         = Def.Ebv;
   Ndof        = Def.Ndof;
elseif (nargin==6),
   Ndof        = Def.Ndof;
elseif (nargin==7),
   % do nothing
else
   error('Illegal number of input arguments');
end


if (isempty(TempVec)),
   TempVec = logspace(3,6,30).';
end

if (isempty(WaveUnits)),
   WaveUnits = Def.WaveUnits;
end
if (isempty(FluxUnits)),
   FluxUnits = Def.FluxUnits;
end

switch lower(FluxUnits)
 case 'f_lambda'
    % do nothing
 otherwise
    error('This flux units is unsupported yet');
end

%--- convert units ---
Spec(:,1) = convert.energy(WaveUnits,'A',Spec(:,1));


%--- Apply extinction ---
Rv = 3.08;
if (Ebv~=0),
   Aw = optical_extinction(Ebv,'B','V',Spec(:,1)./1e4,'C',Rv);
   Spec(:,2) = Spec(:,2).*10.^(0.4.*Aw);
end


Nt = length(TempVec);
Factor = zeros(Nt,1);
RMS    = zeros(Nt,1);
for It=1:1:Nt,
   [~,~,IlA]=black_body(TempVec(It),Spec(:,1),'P');

   %Factor(It) = mean(Spec(:,2)./IlA);
   %[mean(Spec(:,2)./IlA), IlA\Spec(:,2)]
   Factor(It) = IlA\Spec(:,2);
   ErrCL=err_cl(Spec(:,2) - IlA.*Factor(It));
   RMS(It)    = std(Spec(:,2) - IlA.*Factor(It));
   %RMS(It)    = std(log10(Spec(:,2)./IlA.*Factor(It)));
   %RMS(It)    = ErrCL(1,2) - ErrCL(1,1);
end


NoErr = 1;
if (nargout>2),
   [Res.BestRMS,Res.Ind] = min(RMS);

   Extram = find_local_extramum(TempVec,RMS);
   if (isempty(Extram)),
      Res = NaN;
   else
      I      = find(Extram(:,3)>0);
      Extram = Extram(I,:);
      [~,MinI] = min(Extram(:,2));
      Res.BestT = Extram(MinI,1);
      if (isempty(Res.BestT)),
          [MinRMS,MinI] = min(RMS);
          Res.BestT     = TempVec(MinI);
          Res.BestFactor = interp1(TempVec,Factor,Res.BestT,'nearest');
          Res.BestRadius = sqrt(Res.BestFactor.*(Dist.*Pc).^2);
          NoErr = 1;
      else
         Res.BestFactor = interp1(TempVec,Factor,Res.BestT,'cubic');
         Res.BestRadius = sqrt(Res.BestFactor.*(Dist.*Pc).^2);
         NoErr = 0;
      end
   end

   % estimate errors
   % normalize RMS such that \chi^ per dof = 1
   Chi2 = Ndof.*(RMS./min(RMS)).^2;
   %plot(TempVec,Chi2)
   [MinChi2,IndMC] = min(Chi2);
   SubChi2 = Chi2 - MinChi2 - chi2inv(0.68,2);
   if (NoErr == 0),
      [~,IndLE] = min(abs(SubChi2(1:IndMC)));
      [~,IndUE] = min(abs(SubChi2(IndMC:end)));
      Res.TempLowErr  = Res.BestT - TempVec(IndLE);
      Res.TempHighErr = TempVec(IndUE+IndMC) - Res.BestT;
   else
       Res.TempLowErr = Inf;
       Res.TempHighErr = Inf;
   end
else
   Res = [];
end
