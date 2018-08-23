function [Spec,SN_AgeRange]=get_spectra(TemplateName,Ebv,R,Z,Vel,Age,AirMass,Telluric)
%--------------------------------------------------------------------------
% get_spectra function                                           AstroSpec
% Description: Get a template spectrum from spectra library,
%              or alternatively a use supplied spectrum.
%              Optionally, apply redshift, extinction, atmospheric
%              extinction, and Telluric absorptions.
% Input  : - Template name or user supplied spectra [Wave(ang), Intensity]:
%            Templates name can be one of the followings:
%            --- Quasars ---
%            'QSO_LBQS','QSO_FBQS','QSO_FBQS_RL','QSO_FBQS_RQ','QSO_SDSS',
%            'QSO_HBal','QSO_LBal'
%            --- Galaxies ---
%            'Gal_E','Gal_S0','Gal_Sa','Gal_Sb','Gal_Sc','Gal_Sbc'
%            --- StarBurst ---
%            'Gal_StarBurst1','Gal_StarBurst2','Gal_StarBurst3',
%            'Gal_StarBurst4','Gal_StarBurst5','Gal_StarBurst6','Gal_Bulge'
%            --- Stars --- (Silva & Cornel 1992)
%            'O5V','O7B0V','B34V','B6V','A13V','A57V','A8V',
%            'F67V','F89V','G12V','G68V','K4V','K5V','M2V','M4V'
%            'O7B1III','B5III','B9III','A3III','F47III','G04III',
%            'K2III','K4III','K7III','M01III','M3III','M4III','M5III','M6III' 
%            'O8I','B1I','A03I','A79I','F7I','G01I','K35I'
%            --- Stars --- (Pickles 1998)
%            'o5v','o8iii','o9v',
%            'b0i','b0v','b12iii','b1i','b1v','b2ii','b2iv',
%            'b3i','b3iii','b3v','b57v','b5i','b5ii','b5iii','b6iv','b8i','b8v',
%            'b9iii',b9v','ukb0v','ukb12iii','ukb1i','ukb2ii','ukb2iv',
%            'a0i','a0iii','a0iv','a0v','a2i','a2v','a3v',a47iv','a5iii','a5v',
%            'a7iii','a7v','uka0i','uka0iii','uka0iv','uka0v','uka2i','uka2v',
%            'uka3iii','uka3v','uka47iv','uka5v','uka7iii','uka7v',
%            'f02iv','f0i','f0ii','f0iii','f0v','f2ii','f2iii','f2v','f5i',
%            'f5iii','f5iv','f5v','f6v','f8i','f8iv','f8v','rf6v','rf8v',
%            'g0i',g0iii','g0iv','g0v','g2i','g2iv','g2v','g5i','g5ii','g5iii',
%            'g5iv','g5v','g8i','g8iii','g8iv','g8v','rg0v','rg5iii','rg5v',
%            'k01ii','k0iii','k0iv','k0v','k1iii','k1iv','k2i','k2iii','k2v',
%            'k34ii','k3i','k3iii','k3iv','k3v','k4i','k4iii','k4v','k5iii',
%            'k5v','k7v','rk0iii','rk0v','rk1iii','rk2iii','rk3iii','rk4iii',
%            'rk5iii',
%            'm0iii','m0v','m10iii','m1iii','m1v','m2i','m2iii','m2p5v','m2v',
%            'm3ii','m3iii','m3v','m4iii','m4v','m5iii','m5v','m6iii','m6v',
%            'm7iii','m8iii','m9iii'
%            --- SN (time dependent) ---
%            Not available in the public version
%            'SN1984l','SN1987l','SN1990b','SN1990u',
%            'SN1991ar','SN1991bg','SN1991t','SN1992h','SN1993j',
%            'SN1994ak','SN1994d','SN1994i','SN1994y','SN1995d',
%            'SN1996cb','SN1998bp','SN1998bw','SN1998de',
%            'SN1998dt','SN1998es','SN1998s',
%            'SN1999da','SN1999di','SN1999dk','SN1999dn','SN1999ee',
%            'SN1999em','SN2000cx','SN2001x','SN2002ap'
%            'SN1996L','SN1999ex','SN1997cy','SN1998S','SN1991T'
%
%          - Redened spectra using exitinction E(B-V), default is [0 0].
%            The first value is for E(B-V) in the source redshift,
%            while the second value is for the E(B-V) in the observer frame.
%          - Redened spectra using extinction law R, default is [3.1 3.1].
%            The first value is for R in the source redshift, while the second
%            value is for the R in the observer frame.
%          - Redshift of spectra, default is 0.
%          - Relative velocity of spectra, default is 0 -
%            if given override redshift.
%          - SN Age (from maximum light), default is 0.
%          - AirMass extinction to applay to spectrum.
%            Default is 0 (no extinction).
%            Atmospheric extinction used is based on KPNO extinction
%            (source: IRAF).
%            Note this is available only for the 3200A-10400A range.
%          - Add Telluric absorption to spectra,
%            where the value indicates the depth of the 7600A Telluric
%            absorption relative to the template spectra. 0 means nothing to
%            add and 1 mean 71% absorbtion at 7606A, and 34% at 6875A.
%            Default is 0.
% Output : - Spectra [Wavelength(Ang), Intensity].
%          - SN Age range [min, max].
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Spec]=get_spectra('QSO_LBQS');
% Reliable: 2
%--------------------------------------------------------------------------

ColW = 1;
ColI = 2;
InterpMethod = 'linear';

if (nargin==1),
   Ebv      = [0 0];
   R        = [3.1 3.1];
   Z        = 0;
   Vel      = [];
   Age      = 0;
   AirMass  = 0;
   Telluric = 0;
elseif (nargin==2),
   R        = [3.1 3.1];
   Z        = 0;
   Vel      = [];
   Age      = 0;
   AirMass  = 0;
   Telluric = 0;
elseif (nargin==3),
   Z        = 0;
   Vel      = [];
   Age      = 0;
   AirMass  = 0;
   Telluric = 0;
elseif (nargin==4),
   Vel      = [];
   Age      = 0;
   AirMass  = 0;
   Telluric = 0;
elseif (nargin==5),
   Age      = 0;
   AirMass  = 0;
   Telluric = 0;
elseif (nargin==6),
   AirMass  = 0;
   Telluric = 0;
elseif (nargin==7),
   Telluric = 0;
elseif (nargin==8),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (length(Ebv)~=2| length(R)~=2),
   error('E(B-V) and R must be given in both source and observer frame');
end

SN_AgeRange = [NaN NaN];
%--- Read Spectra ---
if (ischar(TemplateName)==1),
   if (strcmp(TemplateName(1:2),'SN')==1),
      %--- SN spectrum in structure ---
      load(sprintf('%s.mat',TemplateName));
      Spec = interp_sn_spec(SN,Age);
      SN_AgeRange = [min(SN.day), max(SN.day)];
   else
      %--- QSO/Gal/Star spectra --- 
      % Spec = load(sprintf('%s.mat',TemplateName));
      Iqso = findstr(TemplateName,'QSO');
      Igal = findstr(TemplateName,'Gal');
      if (isempty(Iqso) && isempty(Igal)),
         if (upper(TemplateName(1))==TemplateName(1)),
            % stars (StellarSpec)
            Spec = load(sprintf('%s.txt',TemplateName));
            Spec(:,ColW) = Spec(:,ColW).*10;   % nanometer to ang
         else
            % Pickeles stars (StellarSpec2)
            Spec = load(sprintf('%s.mat',TemplateName));
            Spec = getfield(Spec,TemplateName);
            % already in Ang.
         end
      else
         % QSO/Gal
         Spec = load(sprintf('%s.txt',TemplateName));
      end
   end 

else
   Spec = TemplateName;
   if (isstruct(Spec)==1),
      % assume SN Spectrum 
      Spec = interp_sn_spec(Spec,Age);
      SN_AgeRange = [min(SN.day), max(SN.day)];
   end
end

%--- convert stellar spectra to Ang ---

%--- Redened spectra in source frame ---
ONES = ones(size(Spec(:,ColW)));
%if (Ebv(1)>0),
   A_w = optical_extinction(Ebv(1),'B','V',Spec(:,ColW)./1e4,'C',R(1).*ONES);    % use  Cardelli et al. model
   % convert magnitude to flux factor
   Spec(:,ColI) = Spec(:,ColI) .* 10.^(-0.4.*A_w); 
%end

%---Redshift spectra ---
if (isempty(Vel)),
   % use redshift
else
   % use velocity
   Z = vel2shift(Vel);
end
if (Z~=0),
    
   Spec = shift_spec(Spec,Z,'w','f');    % applay redshift
end

%--- Redened spectra in observer frame ---
%if (Ebv(2)>0),
   A_w = optical_extinction(Ebv(2),'B','V',Spec(:,ColW)./1e4,'C',R(2).*ONES);    % use  Cardelli et al. model
   % convert magnitude to flux factor
   Spec(:,ColI) = Spec(:,ColI) .* 10.^(-0.4.*A_w); 
%end



%-------------------------------------
%--- Applay atmospheric extinction ---
%-------------------------------------
if (AirMass>0),
   load kpnoextinct.mat;

   AtmosphericExtin = interp1(kpnoextinct(:,1), 2.512.^(-AirMass.*kpnoextinct(:,2)), Spec(:,1), InterpMethod);

   Spec(:,2) = Spec(:,2).*AtmosphericExtin;

end

%-------------------------------
%--- Add Telluric absorption ---
%-------------------------------
if (Telluric>0),
   load Telluric_6875.mat
   load Telluric_7606.mat

   Telluric_6875(:,2) = 1./(1 + Telluric.*(1./Telluric_6875(:,2) - 1));
   Telluric_7606(:,2) = 1./(1 + Telluric.*(1./Telluric_7606(:,2) - 1));

   T6875 = interp1(Telluric_6875(:,1), Telluric_6875(:,2), Spec(:,1), InterpMethod);
   T7606 = interp1(Telluric_7606(:,1), Telluric_7606(:,2), Spec(:,1), InterpMethod);

   I6875 = find(isnan(T6875)==1);
   I7606 = find(isnan(T7606)==1);

   T6875(I6875) = 1;
   T7606(I7606) = 1;

   AllTelluric  = T6875.*T7606;

   Spec(:,2) = Spec(:,2).*AllTelluric;

   if (Telluric>1.4),
      Ineg = find(Spec(:,2)<0);
      Spec(Ineg,2) = 0;
   end

end
   







%-----------------------------------------
function Spec=interp_sn_spec(SN,Age);
%-----------------------------------------

      %--- check if Age is within range ---
      if (Age>=min(SN.day) & Age<=max(SN.day)),
         % range ok - interpolate spectrum
         Nage  = length(SN.day);
         Nw    = length(SN.wl);
         Days  = ones(Nw,1)*SN.day;

         Del = SN.day - Age;
         In = find(Del<0);
         [Max,Ind] = max(Del(In));
         Ind1 = In(Ind);
         Ind2 = Ind1+1; 
         Day1 = SN.day(Ind1);
         Day2 = SN.day(Ind2);
         WeightD1 = abs(Age-Day1)./(Day2-Day1);
         WeightD2 = abs(Age-Day2)./(Day2-Day1);
         Spec  = [SN.wl, WeightD1.*SN.spec(:,Ind1) + WeightD2.*SN.spec(:,Ind2)];

      else
         disp('SN Age is out of range');
         Spec = [SN.wl,SN.spec(:,1)].*NaN;
      end


