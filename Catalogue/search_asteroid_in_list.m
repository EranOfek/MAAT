function Found=search_asteroid_in_list(ImagesList,ColList,FOV_RA,FOV_Dec,AstName,GeodPos);
%----------------------------------------------------------------------------
% search_asteroid_in_list function                                 Catalogue
% Description: Given a list of images (dates and positions) check if a given
%              minor planet may be observed in one of the images.
% Input  : - List of all images.
%          - Index of the following columns in the list of images:
%            [JD, RA, Dec]
%            assuming JD is in days, and coordinates are in degrees.
%          - RA field of view of each image [deg].
%          - Dec field of view of each image [deg].
%          - Asteroid name or number.
%          - Observatory geodetic position [lon(deg), lat(deg), height(km)].
% Output : - Structure containing the matches - with the following fields:
%            .Image - the matched lines from the image catalog.
%            .RA    - Vector of matched RA of the minor planet.
%            .Dec   - Vector of matched Dec of the minor planet.
%            .CooJD - Vector of matched JD for the minor planet position.
%            .Ephem - Vector of structures containing the Horizons
%                     ephemerides for each match.
% Tested : Matlab 7.11
%     By : Eran O. Ofek             December 2010
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: Found=search_asteroid_in_list(PTF_ImagesCatalog,[1 3 4],2048./3600,4096./3600,596,[-116.8599, 33.3574, 2]);
%----------------------------------------------------------------------------

%ImagesList = PTF_ImagesCatalog;
%ColList    = [1 3 4];
%FOV_RA     = 2048./3600;
%FOV_Dec    = 4096./3600;
%AstName    = 596;
%GeodPos    = [-116.8599, 33.3574, 2];

RAD = 180./pi;

Col.JD  = ColList(1);
Col.RA  = ColList(2);
Col.Dec = ColList(3);

% sort catalog by declination
ImagesList = sortrows(ImagesList,Col.Dec);

Nim        = size(ImagesList,1);
% convert FOV from scalar to vector
if (numel(FOV_RA)==1),
   FOV_RA   = FOV_RA.*ones(Nim,1);
end
if (numel(FOV_Dec)==1),
   FOV_Dec   = FOV_Dec.*ones(Nim,1);
end

% range of dates in images list
ExtraDays = 10;
MinJD = floor(min(ImagesList(:,Col.JD)))-ExtraDays;
MaxJD = ceil(max(ImagesList(:,Col.JD)))+ExtraDays+1;


% generate ephemerids for asteroid:
TimeStep = 1;   % [hr]
Ephem = get_horizons(MinJD,MaxJD,AstName,'Geod',GeodPos,'StepSize',TimeStep,'StepUnit','h','Type','p');
Ndays = length(Ephem.J_RA);

DistThresh = max([FOV_RA;FOV_Dec]).*1.05.*sqrt(2)./2./RAD;  % [radians]

Counter = 0;
Found.Image = zeros(0,size(ImagesList,2));
Found.RA    = zeros(0,1);
Found.Dec   = zeros(0,1);
Found.CooJD = zeros(0,1);
Found.Ephem = zeros(0,1);
for Idays=1:1:Ndays-1,
   [Lines,DistL,PAL]=cat_search(ImagesList(:,[Col.RA, Col.Dec])./RAD, [1 2], [Ephem.J_RA(Idays), Ephem.J_Dec(Idays)], DistThresh);

   if (~isempty(Lines)),
      % candidate found - compare dates
      DeltaTime = ImagesList(Lines,Col.JD) - Ephem.JD(Idays);
      Ndt       = length(DeltaTime);
      for Idt=1:1:Ndt,
	 if (abs(DeltaTime(Idt))<TimeStep./24)
            % possible match found
   	    ImageInd = Lines(Idt);
            ImageJD  = ImagesList(ImageInd,Col.JD);
            ImageRA  = ImagesList(ImageInd,Col.RA);
            ImageDec = ImagesList(ImageInd,Col.Dec);

            % generate ephemerids at exact time of image
            EphemIm = get_horizons(ImageJD,ImageJD+1,AstName,'Geod',GeodPos,'Type','p');
            [Offset.RA, Offset.Dec] = sphere_offset(EphemIm.J_RA(1),EphemIm.J_Dec(1),ImageRA./RAD,ImageDec./RAD);
            % check if image contains asteroid
            if (abs(Offset.RA)<FOV_RA(ImageInd)./RAD & abs(Offset.Dec)<FOV_Dec(ImageInd)./RAD),
               % Found asteroid in image
  	       Counter = Counter + 1;
               Found.Image(Counter,:) = ImagesList(ImageInd,:);
               Found.RA(Counter)      = EphemIm.J_RA(1);
               Found.Dec(Counter)     = EphemIm.J_Dec(1);
               Found.CooJD(Counter)   = EphemIm.JD(1);
               Found.Ephem(Counter).EpehmIm   = EphemIm;
jd2date(Found.CooJD(Counter))

               % Asteroid found in image - therefore remove image from ImagesList
               ImagesList = Util.array.delete_ind(ImagesList,ImageInd,1);

            end
         end
      end
   end
end
