function RandCoo=cel_coo_rnd(N,varargin);
%------------------------------------------------------------------------------
% cel_coo_rnd function                                               AstroStat
% Description: Generate random coordinates on the celestial sphere.
%              The program can applay matrix rotation, reject
%              coordinates, and generate non uniform coordinates.
%              [preliminary version].
% Input  : - If scalar then: Number of random element to generate,
%            Else: a list of coordinates [Long, Lat, Err, ErrType] (in radians)
%            for which to add random error to each coordinates,
%            using the error column.
%            ErrType is: 0 - for Gaussian errors
%                        1 - for Fisher errors ('Gaussian on a sphere').
%          * Aribtrary number of pairs of arguments:
%            ...,keyword,value,...
%            --- Note the program execute the commands by their order
%                from first to last ---
%            Possible keywords are:
%            'Rotate'       - Applay a rotation matrix to coordinates.
%            'RejLatRange'  - Reject Galactic latitude range [MinLat, MaxLat],
%                             default is [] (don't rejuct).
%            'RejLongRange' - Reject Galactic longitute range [MinLat, MaxLat],
%                             default is [] (don't rejuct).
%                             default is [] (don't rejuct).
%            'RejPoly'      - Cell array of polygons to rejuct.
%                             default is [] (don't rejuct).
%            'ProbLong'     - Galactic longitude probability [Long, Prob],
%                             default is, [], flat probability (i.e., Prob=1).
%            'ProbLat'      - Galactic latitude probability [Lat, Prob],
%                             default is, [], flat probability (i.e., Prob=1).
% Output : - Random coordinates [Long, Lat] in radians.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                  November 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: RandCoo=cel_coo_rnd(1000,'RejLatRange',[-20 20; 85 86]./RAD);
%          RandCoo=cel_coo_rnd(1000,'RejLatRange',[-20 20; 85 86]./RAD,...
%                                   'Rotate',RotMat,'ProbLat',ProbLat);
%          % (1) reject in Galactic
%          % (2) rotate to equatorial
%          % (3) applay probability map
%          RandCoo=cel_coo_rnd(1000,'RejLatRange',[-20 20; 85 86]./RAD,..
%                                   'Rotate',rotm_coo('G'),...
%                                   'ProbLat',[-pi./2 1; +pi./2 0.5]);
% Reliable: 2
%------------------------------------------------------------------------------

RAD          = 180./pi;
Nfactor      = 200;
ColLong      = 1;
ColLat       = 2;
InterpMethod = 'linear';
Col_Err      = 3;
Col_ErrType  = 4;


if (prod(size(N))==1)
   Method = 'Random';
else
   Method = 'List';
end

switch Method
 case 'Random'
    %--- select random coordinates ---
    SkyEffFactor = Nfactor.*(8./(4.*pi./3));
    
    EffectiveN = ceil(SkyEffFactor.*N);
    RandCooDir = rand(EffectiveN,3).*2 - 1;   % Random cosine direction
    % select vector within unit sphere
    I          = find( sum(RandCooDir.^2,2) <= 1 );
    RandCooDir = RandCooDir(I,:);
    % Normalize to unit vector
    RandCooDir = RandCooDir./(sqrt(sum(RandCooDir.^2,2))*ones(1,3));
    RandCoo    = celestial.coo.cosined(RandCooDir);              % random coordinates [radians]
    % convert to 0..2pi range
    I          = find(RandCoo(:,ColLong)<0);
    RandCoo(I,ColLong) = 2.*pi + RandCoo(I,ColLong);

 case 'List'
    %--- Select random coordinates by adding random errors to list ---
    InList = N;
    N      = size(InList,1);



    if (length(find(InList(:,Col_ErrType)==0)).*length(find(InList(:,Col_ErrType)==1))~=0)
       error('In this version: ErrorType must be 1 | 0 but not mixed');
    end

    switch InList(1,Col_ErrType)
     case 0
        % Gaussian errors
        Range   = abs(randn(N,1).*InList(:,Col_Err));
        Azimuth = rand(N,1).*2.*pi;
        [LatOut,LongOut] = reckon(InList(:,ColLat), InList(:,ColLong), Range, Azimuth,'radians');

        RandCoo = [LongOut, LatOut];
     case 1
        % Fisher errors
	error('        not available yet');
     otherwise
        error('Unknown ErrFunCode Option');
    end

 otherwise
    error('Unknown Method Option');
end




%--- Applay rejuction and rotation by their order ---
Narg = length(varargin);
for Iarg=1:2:Narg-1
   switch varargin{Iarg}
    case 'RejLatRange'
       %--- Reject by latitude ---
       RejLatRange = varargin{Iarg+1};
       Nregions = size(RejLatRange,1);
       for I=1:1:Nregions
          Iselected = find(RandCoo(:,ColLat)<RejLatRange(I,1) | RandCoo(:,ColLat)>RejLatRange(I,2));
          RandCoo   = RandCoo(Iselected,:);
       end

    case 'RejLongRange'
       %--- Reject by longitude ---
       RejLongRange = varargin{Iarg+1};
       Nregions = size(RejLongRange,1);
       for I=1:1:Nregions
          Iselected = find(RandCoo(:,ColLong)<RejLongRange(I,1) | RandCoo(:,ColLong)>RejLongRange(I,2));
          RandCoo   = RandCoo(Iselected,:);
       end

    case 'ProbLat'
       %--- Reject by probability map in latitude ---
       ProbLat = varargin{Iarg+1};
       % interpolate probability rejection map
       NewN      = size(RandCoo,1);
       Iselected = find( interp1(ProbLat(:,1),ProbLat(:,2),RandCoo(:,ColLat),InterpMethod) > rand(NewN,1));
       RandCoo   = RandCoo(Iselected,:);

    case 'ProbLong'
       %--- Reject by probability map in longitude ---
       ProbLong = varargin{Iarg+1};
       % interpolate probability rejection map
       NewN      = size(RandCoo,1);
       Iselected = find( interp1(ProbLong(:,1),ProbLong(:,2),RandCoo(:,ColLong),InterpMethod) > rand(NewN,1));
       RandCoo   = RandCoo(Iselected,:);

    case 'RejPoly'
       %--- Reject inside polygons ---
       RejPoly = varargin{Iarg+1};
       Nregions = length(RejPoly);
       for I=1:1:Nregions
          Iselected = find(inspherepolygon(RandCoo(:,ColLong),RandCoo(:,ColLat),RejPoly{I}(:,ColLong),RejPoly{I}(:,ColLat))==0);
          RandCoo   = RandCoo(Iselected,:);
       end

    case 'Rotate'

       %--- Applay rotation matrix ---
       Rotate  = varargin{Iarg+1};
       RandCoo = celestial.coo.cosined( (Rotate*[celestial.coo.cosined(RandCoo).']).');

       I = find(RandCoo(:,ColLong)<0);
       RandCoo(I,ColLong) = 2.*pi + RandCoo(I,ColLong);

    otherwise
      error('Unknown Keyword Option');
   end
end



switch Method
 case 'Random'
    if (size(RandCoo,1)<N)
       error('Number of random coordinates samller then needed - increase SkyEffFactor');
    else
       RandCoo = RandCoo(1:N,:);
    end
 case 'List'
     % do nothing
 otherwise
    error('Unknown Method Option');
end


