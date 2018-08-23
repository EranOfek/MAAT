function [BestOutMag,MinRMS,BestFitTemplate,RMS]=interp_mag(MeasuredMag,InFilters,InABV,OutFilters,OutABV,Z,Ebv,TemplateType,AirMass,Telluric);
%--------------------------------------------------------------------------
% interp_mag function                                            AstroSpec
% Description: Given magnitude of an object find the best fit spectra
%              (from a library of templates), and calculate the magnitude of
%              the spectra in additional bands.
% Input  : - Input measured magnitudes, and optional errors [Mag, Err].
%            If single column is given, then give equal weight
%            for all measurments.
%          - Cell array of filters normalized transmission,
%            corresponding to the input magnitudes.
%          - Cell array of magnitude types for input mag 'AB' | 'Vega'.
%          - Cell array of filters normalized transmission,
%            in which to calculate the output magnitudes.
%          - Cell array of magnitude types for input mag 'AB' | 'Vega'.
%          - Vector of redshifts to test, default is 0.
%          - Vector of E_{B-V} extinctions to test, default is 0.
%          - Types of templates to use {'Star' | 'StarMS' | 'QSO' | 'Gal'},
%            default is 'Star'.
%          - Two elements vector containing the airmass of atmospheric
%            extinction to applay to a spectrum while
%            calculating the input magnitudes (first element), and
%            while calculating the output magnitudes (second element).
%            Default is [0 0] - no extinction.
%          - Two elements vector containing the strength of the
%            Telluric absorptions (see get_spectra.m)
%            to applay to a spectrum while calculating the input
%            magnitudes (first element), and while calculating
%            the output magnitudes (second element).
%            Default is [0 0] - no Telluric absorption.
% Output : - Vector of output magnitudes.
%          - Best RMS [mag].
%          - Best fit template name.
%          - RMS matrix for all templates.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Needed: Spectra library
% Example: % Fit template and interpolate magnitudes
%          InABV  = {'AB','AB'}; OutABV = {'AB','AB','AB'};
%          MeasuredMag=[19.69;19.24];
%          F475W=get_filter('HST-ACS','F475W');
%          F814W=get_filter('HST-ACS','F814W');
%          g=get_filter('SDSS','g');
%          r=get_filter('SDSS','r');
%          i=get_filter('SDSS','i');
%          [BestOutMag,MinRMS,BestFit,RMS]=interp_mag(MeasuredMag,{F475W.nT{1},F814W.nT{1}},InABV,{g.nT{1},r.nT{1},i.nT{1}},OutABV,0,0)
% Reliable: 1
%--------------------------------------------------------------------------
StatType = 'wmean';    % {'mean' | 'wmean'} - use measured mag weights
R        = 3.08;       % R = A_{V} / E_{B-V}

if (nargin==5),
   Z            = 0;
   Ebv          = 0;
   TemplateType = 'Star';
   AirMass      = [0 0];
   Telluric     = [0 0];
elseif (nargin==6),
   Ebv          = 0;
   TemplateType = 'Star';
   AirMass      = [0 0];
   Telluric     = [0 0];
elseif (nargin==7),
   TemplateType = 'Star';
   AirMass      = [0 0];
   Telluric     = [0 0];
elseif (nargin==8),
   AirMass      = [0 0];
   Telluric     = [0 0];
elseif (nargin==9),
   Telluric     = [0 0];
elseif (nargin==10),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (iscell(InFilters)==0);
   InBands{1} = InFilters;
else
   InBands    = InFilters;
end

if (iscell(OutFilters)==0);
   OutBands{1} = OutFilters;
else
   OutBands    = OutFilters;
end

if (size(MeasuredMag,2)==1),
   MeasuredMag = [MeasuredMag, 0.01.*ones(size(MeasuredMag))];
end


%StarNames = {'a0iii';'a0i';'a0iv';'a0v';'a2i';'a2v';'a3iii';'a3v';'a47iv';'a5iii';'a5v';'a7iii';'a7v';'b0i';'b0v';'b12iii';'b1i';'b1v';'b2ii';'b2iv';'b3iii';'b3i';'b3v';'b57v';'b5iii';'b5ii';'b5i';'b6iv';'b8i';'b8v';'b9iii';'b9v';'f02iv';'f0iii';'f0ii';'f0i';'f0v';'f2iii';'f2ii';'f2v';'f5iii';'f5i';'f5iv';'f5v';'f6v';'f8i';'f8iv';'f8v';'g0iii';'g0i';'g0iv';'g0v';'g2i';'g2iv';'g2v';'g5iii';'g5ii';'g5i';'g5iv';'g5v';'g8iii';'g8i';'g8iv';'g8v';'k01ii';'k0iii';'k0iv';'k0v';'k1iii';'k1iv';'k2iii';'k2i';'k2v';'k34ii';'k3iii';'k3i';'k3iv';'k3v';'k4iii';'k4i';'k4v';'k5iii';'k5v';'k7v';'m0iii';'m0v';'m10iii';'m1iii';'m1v';'m2iii';'m2i';'m2p5v';'m2v';'m3iii';'m3ii';'m3v';'m4iii';'m4v';'m5iii';'m5v';'m6iii';'m6v';'m7iii';'m8iii';'m9iii';'o5v';'o8iii';'o9v';'rf6v';'rf8v';'rg0v';'rg5iii';'rg5v';'rk0iii';'rk0v';'rk1iii';'rk2iii';'rk3iii';'rk4iii';'rk5iii';'uka0iii';'uka0i';'uka0iv';'uka0v';'uka2i';'uka2v';'uka3iii';'uka3v';'uka47iv';'uka5iii';'uka5v';'uka7iii';'uka7v';'ukb0i';'ukb0v';'ukb12iii';'ukb1i';'ukb1v';'ukb2ii';'ukb2iv'};


StarNames={'uka0iii';'uka0i';'uka0iv';'uka0v';'uka2i';'uka2v';'uka3iii';'uka3v';'uka47iv';'uka5iii';'uka5v';'uka7iii';'uka7v';'ukb0i';'ukb0v';'ukb12iii';'ukb1i';'ukb1v';'ukb2ii';'ukb2iv';'ukb3iii';'ukb3i';'ukb3v';'ukb57v';'ukb5iii';'ukb5ii';'ukb5i';'ukb6iv';'ukb8i';'ukb8v';'ukb9iii';'ukb9v';'ukf02iv';'ukf0iii';'ukf0ii';'ukf0i';'ukf0v';'ukf2iii';'ukf2ii';'ukf2v';'ukf5iii';'ukf5i';'ukf5iv';'ukf5v';'ukf6v';'ukf8i';'ukf8iv';'ukf8v';'ukg0iii';'ukg0i';'ukg0iv';'ukg0v';'ukg2i';'ukg2iv';'ukg2v';'ukg5iii';'ukg5ii';'ukg5i';'ukg5iv';'ukg5v';'ukg8iii';'ukg8i';'ukg8iv';'ukg8v';'ukk01ii';'ukk0iii';'ukk0iv';'ukk0v';'ukk1iii';'ukk1iv';'ukk2iii';'ukk2i';'ukk2v';'ukk34ii';'ukk3iii';'ukk3i';'ukk3iv';'ukk3v';'ukk4iii';'ukk4i';'ukk4v';'ukk5iii';'ukk5v';'ukk7v';'ukm0iii';'ukm0v';'ukm10iii';'ukm1iii';'ukm1v';'ukm2iii';'ukm2i';'ukm2p5v';'ukm2v';'ukm3iii';'ukm3ii';'ukm3v';'ukm4iii';'ukm4v';'ukm5iii';'ukm5v';'ukm6iii';'ukm6v';'ukm7iii';'ukm8iii';'ukm9iii';'uko5v';'uko8iii';'uko9v';'ukrf6v';'ukrf8v';'ukrg0v';'ukrg5iii';'ukrg5v';'ukrk0iii';'ukrk0v';'ukrk1iii';'ukrk2iii';'ukrk3iii';'ukrk4iii';'ukrk5iii';'ukwf5v';'ukwf8v';'ukwg0v';'ukwg5iii';'ukwg5v';'ukwg8iii';'ukwk0iii';'ukwk1iii';'ukwk2iii';'ukwk3iii';'ukwk4iii'};

StarNamesMS={'uka0v';'uka2v';'uka3v';'uka5v';'uka7v';'ukb0v';'ukb1v';'ukb3v';'ukb57v';'ukb8v';'ukb9v';'ukf0v';'ukf2v';'ukf5v';'ukf6v';'ukf8v';'ukg0v';'ukg2v';'ukg5v';'ukg8v';'ukk0v';'ukk2v';'ukk3v';'ukk4v';'ukk5v';'ukk7v';'ukm0v';'ukm1v';'ukm2p5v';'ukm2v';'ukm3v';'ukm4v';'ukm5v';'ukm6v';'uko5v';'uko9v';'ukrf6v';'ukrf8v';'ukrg0v';'ukrg5v';'ukrk0v';'ukwf5v';'ukwf8v';'ukwg0v';'ukwg5v'};





QSONames  = {'QSO_LBQS';'QSO_FBQS';'QSO_FBQS_RL';'QSO_FBQS_RQ';'QSO_HBal';'QSO_LBal'};

GalNames  = {'Gal_E';'Gal_S0';'Gal_Sa';'Gal_Sb';'Gal_Sc';'Gal_Sbc';'Gal_StarBurst1';'Gal_StarBurst2';'Gal_StarBurst3';'Gal_StarBurst4';'Gal_StarBurst5';'Gal_StarBurst6';'Gal_Bulge'};


switch TemplateType
 case 'Star'
    ListNames = StarNames;
 case 'StarMS'
    ListNames = StarNamesMS;
 case 'QSO'
    ListNames = QSONames;
 case 'Gal'
    ListNames = GalNames;
 otherwise
    error('Unknown TemplateType Option');
end

Nobj      = length(ListNames);
Nz        = length(Z);
Nebv      = length(Ebv);
RMS       = zeros(Nobj,Nz,Nebv);
Mean      = zeros(Nobj,Nz,Nebv);
WMean     = zeros(Nobj,Nz,Nebv);
InMagMat  = zeros(Nobj,Nz,Nebv,length(InBands));
OutMagMat = zeros(Nobj,Nz,Nebv,length(OutBands));
for Iobj=1:1:Nobj,
   for Iz=1:1:Nz,
      for Iebv=1:1:Nebv,

         [Spec] = get_spectra(ListNames{Iobj},[0 Ebv(Iebv)],[R R],Z(Iz),[],0,AirMass(1),Telluric(1));

         Diff = zeros(length(InBands),1);
         Err  = zeros(length(InBands),1);
         for Iband=1:1:length(InBands),
            %[InMag,Dmag,Extrap_frac] = basic_synthetic_photometry(Spec,InBands{Iband},InABV{Iband},[],[0 1]);
            [InMag,Extrap_frac] = synphot(Spec,InBands{Iband},[],InABV{Iband});
            InMagMat(Iobj,Iz,Iebv,Iband) = InMag;
            Diff(Iband) = MeasuredMag(Iband,1) - InMag;
            Err(Iband) = MeasuredMag(Iband,2);   % Measurment error
            if (Extrap_frac>0.01),
	       Diff(Iband) = NaN;
            end
         end

         RMS(Iobj,Iz,Iebv)   = std(Diff);
         Mean(Iobj,Iz,Iebv)  = mean(Diff);
         WMean(Iobj,Iz,Iebv) = wmean([Diff,Err]);

         if (AirMass(1)==AirMass(2) && Telluric(1)==Telluric(2)),
            % do nothing - Use Spec as is
	 else
            % calculate new Spec
            [Spec] = get_spectra(ListNames{Iobj},[0 Ebv(Iebv)],[R R],Z(Iz),[],0,AirMass(2),Telluric(2));
         end

         for Iband=1:1:length(OutBands),
            %[OutMag,Dmag,Extrap_frac] = basic_synthetic_photometry(Spec,OutBands{Iband},OutABV{Iband},[],[0 1]);
            [OutMag,Extrap_frac] = synphot(Spec,OutBands{Iband},[],OutABV{Iband});

            OutMagMat(Iobj,Iz,Iebv,Iband)       = OutMag;
            OutExtrapMagMat(Iobj,Iz,Iebv,Iband) = Extrap_frac;
         end
       
      end
   end
end


[MinRMS,MinInd] = minnd(RMS);
if (length(MinInd)==2),
   switch StatType
    case 'mean'
       BestOutMag = OutMagMat(MinInd(1),MinInd(2),:) + Mean(MinInd(1),MinInd(2));
    case 'wmean'
       BestOutMag = OutMagMat(MinInd(1),MinInd(2),:) + WMean(MinInd(1),MinInd(2));
    otherwise
       error('Unknown StatType Option');
   end
else
   switch StatType
    case 'mean'
       BestOutMag = OutMagMat(MinInd(1),MinInd(2),MinInd(3),:) + Mean(MinInd(1),MinInd(2));
    case 'wmean'
       BestOutMag = OutMagMat(MinInd(1),MinInd(2),MinInd(3),:) + WMean(MinInd(1),MinInd(2));
    otherwise
       error('Unknown StatType Option');
   end
end

BestOutMag = squeeze(BestOutMag);

BestFitTemplate = ListNames(MinInd(1));

%I = find(RMS<0.07)
%OutMagMat(I,MinInd(2),:)
%Mean(I,MinInd(2))


