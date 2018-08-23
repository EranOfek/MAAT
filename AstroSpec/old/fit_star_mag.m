function Best=fit_star_mag(Filters,Mags)

% Input  : - Three column cell array of filters in which the magnitudes are specified.
%            The columns are: family name, filter name, magnitude system.
%            e.g., {'SDSS','g','AB'; '2MASS','J','Vega'; 'POSS-II','Red','Vega'}
%          - Matrix of magnitudes in the various filters.
%            Each column corresponds to one object and each row to a filter
%            (corresponding to the filter in the Filters cell array).

%          - Optional vector of E_{B-V} to fit... not implemented

RAD = 180./pi;




StarNames={'uka0iii';'uka0i';'uka0iv';'uka0v';'uka2i';'uka2v';'uka3iii';'uka3v';'uka47iv';'uka5iii';'uka5v';'uka7iii';'uka7v';'ukb0i';'ukb0v';'ukb12iii';'ukb1i';'ukb1v';'ukb2ii';'ukb2iv';'ukb3iii';'ukb3i';'ukb3v';'ukb57v';'ukb5iii';'ukb5ii';'ukb5i';'ukb6iv';'ukb8i';'ukb8v';'ukb9iii';'ukb9v';'ukf02iv';'ukf0iii';'ukf0ii';'ukf0i';'ukf0v';'ukf2iii';'ukf2ii';'ukf2v';'ukf5iii';'ukf5i';'ukf5iv';'ukf5v';'ukf6v';'ukf8i';'ukf8iv';'ukf8v';'ukg0iii';'ukg0i';'ukg0iv';'ukg0v';'ukg2i';'ukg2iv';'ukg2v';'ukg5iii';'ukg5ii';'ukg5i';'ukg5iv';'ukg5v';'ukg8iii';'ukg8i';'ukg8iv';'ukg8v';'ukk01ii';'ukk0iii';'ukk0iv';'ukk0v';'ukk1iii';'ukk1iv';'ukk2iii';'ukk2i';'ukk2v';'ukk34ii';'ukk3iii';'ukk3i';'ukk3iv';'ukk3v';'ukk4iii';'ukk4i';'ukk4v';'ukk5iii';'ukk5v';'ukk7v';'ukm0iii';'ukm0v';'ukm10iii';'ukm1iii';'ukm1v';'ukm2iii';'ukm2i';'ukm2p5v';'ukm2v';'ukm3iii';'ukm3ii';'ukm3v';'ukm4iii';'ukm4v';'ukm5iii';'ukm5v';'ukm6iii';'ukm6v';'ukm7iii';'ukm8iii';'ukm9iii';'uko5v';'uko8iii';'uko9v';'ukrf6v';'ukrf8v';'ukrg0v';'ukrg5iii';'ukrg5v';'ukrk0iii';'ukrk0v';'ukrk1iii';'ukrk2iii';'ukrk3iii';'ukrk4iii';'ukrk5iii';'ukwf5v';'ukwf8v';'ukwg0v';'ukwg5iii';'ukwg5v';'ukwg8iii';'ukwk0iii';'ukwk1iii';'ukwk2iii';'ukwk3iii';'ukwk4iii'};


StarNames={'uka0v';'uka2v';'uka3v';'uka5v';'uka7v';'ukb0v';'ukb1v';'ukb3v';'ukb8v';'ukb9v';'ukf0v';'ukf2v';'ukf5v';'ukf6v';'ukf8v';'ukg0v';'ukg2v';'ukg5v';'ukg8v';'ukk0v';'ukk2v';'ukk3v';'ukk4v';'ukk5v';'ukk7v';'ukm0v';'ukm1v';'ukm2v';'ukm3v';'ukm4v';'ukm5v';'ukm6v';'uko5v';'uko9v'};



FID = fopen('synphot_short','r');
StarInfo = textscan(FID,'%f %f %f %f %f %f %s\n','CommentStyle','%');
fclose(FID);
Col.Teff = 1;
Col.Metal = 2;
Col.Mbol = 3;
Col.BCv = 4;
Col.BCi = 5;
Col.BCk = 6;
Col.Sp = 7;


Nfilt  = size(Filters,1);
Ntempl = length(StarNames);
Nmags  = size(Mags,2);


RMS = zeros(Ntempl,Nmags);
for Itempl=1:1:Ntempl,
   S = load(sprintf('%s.mat',StarNames{Itempl}));
   FN = fieldnames(S);
   Spec = getfield(S,FN{1});
   Data(Itempl).Name = StarNames{Itempl};
   Data(Itempl).Spec = Spec;

   Data(Itempl).SynMag = zeros(Nfilt,1).*NaN;
   for Ifilt=1:1:Nfilt,
      Data(Itempl).SynMag(Ifilt) = synphot(Spec,Filters{Ifilt,1},Filters{Ifilt,2},Filters{Ifilt,3});
   end

   Data(Itempl).SynMagV = synphot(Spec,'Johnson','V','Vega');
   % get additional information on star
   I=find(strcmp(StarInfo{Col.Sp},StarNames{Itempl}(3:end))==1);
   Data(Itempl).Teff   = 10.^StarInfo{Col.Teff}(I);
   Data(Itempl).metal  = StarInfo{Col.Metal}(I);
   Data(Itempl).Mbol   = StarInfo{Col.Mbol}(I);
   Data(Itempl).BCv    = StarInfo{Col.BCv}(I);
   Data(Itempl).BCi    = StarInfo{Col.BCi}(I);
   Data(Itempl).BCk    = StarInfo{Col.BCk}(I);
   Data(Itempl).ModelV = Data(Itempl).Mbol - Data(Itempl).BCv;  % Model V-band mag.

   % fit star
   RMS(Itempl,:) = std(bsxfun(@minus,Mags,Data(Itempl).SynMag));

end

% select best fits
[MinRMS,MinInd] = min(RMS);
Best.Teff = [Data(MinInd).Teff];
Best.Mbol = [Data(MinInd).Mbol];
Best.ModelV = [Data(MinInd).ModelV];
Best.SynMagV = [Data(MinInd).SynMagV];
Best.Name    = [Data(MinInd).Name];
Best.RMS     = min(RMS);

BestSynMag = [Data(MinInd).SynMag] - (Best.SynMagV - Best.ModelV)

Best.Dist = 10.*mean(10.^(0.2.*(Mags-BestSynMag)));

SolLum = 3.839e33;
Sigma  = get_constant('sigma');
Pc     = get_constant('pc');
Best.AngRad = sqrt( SolLum.*10.^(0.4.*(4.74-Best.Mbol))./(4.*pi.*Sigma.*(Best.Dist.*Pc).^2.*Best.Teff.^4 ));
Best.AngRad = Best.AngRad.*RAD.*3600;

