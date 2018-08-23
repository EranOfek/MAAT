function [Cube,SegmentsInfo]=sedm_segments2cube(varargin)
%--------------------------------------------------------------------------
% sedm_segments2cube function                                         SEDM
% Description: Given a SegmentsInfo structure that contains in each
%              element the wavelength calibrated spectrum of each
%              spexcell, generate a cube of flux values as a function
%              of X position, Y position and wavelength.
% Input  : * Arbitrary number of pairs of ...,key,val,.. arguments.
%            The following keywords are available:
%            'SI'    - Input SegmentsInfo structure array as returned by
%                      sedm_wavecalib.m.
%                      If a string is provided then will assume this is
%                      a mat file containing a structure named SegmentsInfo.
%                      This parameter must be provided.
%            'OutFile'- An optional output FITS file name containing
%                      the cube. If empty, then will not write a FITS
%                      file. Default is empty.
%            'Map'   - A structure array containing the mapping between
%                      the index of the SegmentsInfo elemnts and their
%                      spexcell X, Y coordinates (i.e., on sky image
%                      coordinates). This structure should contain the
%                      following fields:
%                      .X   - Spexcell X coordinate [arcsec].
%                      .Y   - Spexcell Y coordinate [arcsec].
%                      .Ind - Index of SegmentsInfo corresponding to
%                             the X/Y coordinates.
%                      Alternatively this can be a string containing
%                      a mat file which contains the Map structure.
%                      Default is 'sedm_MAP.mat'.
%            'Wave'  - A vector of wavelength [nm units] at which
%                      the spexcell spectra will be resampled.
%                      Default is logspace(log10(320),log10(1000),200).'.
%            'Filter'- Run a filter on 1D spectrum before resampling.
%                      This is mainly for increasing S/N.
%                      The following options are available:
%                      'none'     - No filtering (default).
%                      'medfilt1' - run medfilt1 (median filter) on
%                                   spectrum. Parameters are supplied
%                                   using the FilterPars keyword.
%                      'sgolay'   - run sgolayfilt.m
%                                   (Savitzky-Golay filter).                                   
%                      'runmean1'  - run runmean1.m (runninmg mean)
%            'FilterPars' - A cell array of additional parameters to
%                      pass to filters. Default is empty cell {}.
%                      See medfilt1.m or sgolayfilt.m for options.
%                      For run
%            'SpecField'- String containing the name of the spectrum
%                      field in the SegmentsInfo structure array.
%                      Default is 'SpexSpecCenter'.
%            'WaveField'- String containing the name of the wavelength
%                      field in the SegmentsInfo structure array.
%                      Default is 'BestWaveCalib'.
%            'Interp1Method' - Interpolation method for 1D resampling.
%                      See interp1.m for options. Default is 'cubic'.
%            'Interp2Method' - Interpolation method for 2D resampling.
%                      See interp2.m for options. Default is 'cubic'.
%            'CubeInterp' - Cube resampling method. The following options:
%                      'griddata' - use griddata.m function (default).
%                      'tri' -Use the Delaunay triangulation.
%            'VecOnSkyX' - Resampling grid in X direction [arcsec].
%                      Default is (-13:0.5:17.5).
%            'VecOnSkyY' - Resampling grid in Y direction [arcsec].
%                      Default is (-6:0.5:30.5).
%            'GoodRangeX' - Use spexcell within the following X range.
%                      Default is [150 1900].
%            'GoodRangeY' - Use spexcell within the following Y range.
%                      Default is [150 1950].
%            'SpexImSize' - Spexcells image size [X,Y]. Default is [60 60].
% Output : - A structure containing the cube. The following fields are
%            available.
%            .Cube - Cube of [Wavelength (nm), X (pix), Y (pix)].
%            .Wave - Calibrated wavelength [nm].
%            .Xas  - Relative X position [arcsec].
%            .Yas  - Relative X position [arcsec].
%          - The SegmentsInfo structure array with the following new
%            fields:
%            .ResampSpec - The spectrum resampled at the requested
%                          wavelength [Wave, Flux].
%            .OnSkyX     - Spexcell on sky X position [arcsec].
%            .OnSkyY     - Spexcell on sky Y position [arcsec].
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Aug 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Cube,SegmentsInfo]=sedm_segments2cube('SI',SegmentsInfo);
% Reliable: 2
%--------------------------------------------------------------------------

DefV.SI                        = [];
DefV.OutFile                   = [];
DefV.Map                       = 'sedm_MAP.mat';
DefV.Wave                      = logspace(log10(360),log10(950),200).';
DefV.SpecField                 = 'SpexSpecFit';
DefV.WaveField                 = 'WaveCalib';
DefV.Interp1Method             = 'cubic';

DefV.CubeInterp                = 'griddata';   % {'griddata'|'tri'} The Delaunay triangulation
DefV.Interp2Method             = 'cubic';

DefV.Filter                    = 'none'; %medfilt1';
DefV.FilterPars                = {3};
DefV.SpexImSize                = [60 60];
DefV.VecOnSkyX                 = (-13:0.5:17.5);  %(-18:0.5:18);
DefV.VecOnSkyY                 = (-6:0.5:30.5);   %(-13:0.5:13);
DefV.GoodRangeX                = [150 1900];
DefV.GoodRangeY                = [100 1950];

%DefV.MultiThred                = 1;
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

%FunDispLam = @(Pix,Scale,FunPar0,FunPar1) InPar.FunPar0.*(Scale .* InPar.FunPar1).^Pix;

if (ischar(InPar.SI)),
    SegmentsInfo = load2(InPar.SI);
else
   SegmentsInfo = InPar.SI;
   InPar = rmfield(InPar,'SI');
end

if (isempty(InPar.Map)),
    % generate Map
    error('Map generation doesnt exist yet');
end
if (ischar(InPar.Map)),
    InPar.Map = load2(InPar.Map);
end


% % to biuld Map
% clear Map
% load MeanPos_eran.txt
% load MeanPos_26082014.txt
% MeanPos_eran = MeanPos_26082014;
% [~,Iu] = unique(MeanPos_eran(:,1));
% MeanPos_eran = MeanPos_eran(Iu,:);
% % X
% Cell = num2cell(MeanPos_eran(:,4));
% N = size(MeanPos_eran,1);
% [Map(1:1:N).X]=deal(Cell{:});
% % Y
% Cell = num2cell(MeanPos_eran(:,5));
% [Map(1:1:N).Y]=deal(Cell{:});
% % Ind
% Cell = num2cell(MeanPos_eran(:,1));
% [Map(1:1:N).Ind]=deal(Cell{:});
% [~,I] = unique([Map.Ind]);
% [Map1(1:1:length(I)).X] = deal(Map(I).X);
% [Map1(1:1:length(I)).Y] = deal(Map(I).Y);
% [Map1(1:1:length(I)).Ind] = deal(Map(I).Ind);
% save sedm_MAP.mat Map


% go over all Spexcell and generate equaly sampled (wavelength) vectors
Nseg  = length(SegmentsInfo);
CellNaN = num2cell(ones(1,Nseg).*NaN);
All = struct('X',CellNaN,'Y',CellNaN,'Wave',CellNaN,'Int',CellNaN);

Nwave = length(InPar.Wave);
%Cube  = zeros(Nwave,InPar.SpexImSize(1),InPar.SpexImSize(2)).*NaN;
for Iseg=1:1:Nseg,
    SegmentsInfo(Iseg).OnSkyX = NaN;
    SegmentsInfo(Iseg).OnSkyY = NaN;
    
    if (SegmentsInfo(Iseg).MeanX>InPar.GoodRangeX(1) && ...
        SegmentsInfo(Iseg).MeanX<InPar.GoodRangeX(2) && ...
        SegmentsInfo(Iseg).MeanY>InPar.GoodRangeY(1) && ...
        SegmentsInfo(Iseg).MeanY<InPar.GoodRangeY(2)),
        % run Filter
        switch lower(InPar.Filter)
            case 'none'
                % do nothing
                FiltSpec = SegmentsInfo(Iseg).(InPar.SpecField);
            case 'medfilt1'
                FiltSpec = medfilt1(SegmentsInfo(Iseg).(InPar.SpecField),InPar.FilterPars{:});
            case 'sgolay'
                FiltSpec = sgolayfilt(SegmentsInfo(Iseg).(InPar.SpecField),InPar.FilterPars{:});
            case 'runmean1'
                FiltSpec = runmean1(SegmentsInfo(Iseg).(InPar.SpecField),InPar.FilterPars{:});
            otherwise
                error('Unknown Filter option');
        end

        if (isempty(SegmentsInfo(Iseg).(InPar.WaveField))),
            % no wavelength solution
            SegmentsInfo(Iseg).ResampSpec = [];
            %SegmentsInfo(Iseg).ResampWave = [];
            MapInd = [InPar.Map.Ind]==Iseg;
            SegmentsInfo(Iseg).Map        = InPar.Map(MapInd);
        else
           W = SegmentsInfo(Iseg).(InPar.WaveField);
           Inn = ~isnan(W) & ~isnan(FiltSpec);

           if (isempty(find(Inn,1))),
               ResampSpec = NaN.*ones(size(InPar.Wave));
           else
              ResampSpec = interp1(W(Inn),...
                                   FiltSpec(Inn),...
                                   InPar.Wave,...
                                   InPar.Interp1Method);
           end

           SegmentsInfo(Iseg).ResampSpec = [InPar.Wave, ResampSpec]; 
           %SegmentsInfo(Iseg).ResampWave = InPar.Wave;



           MapInd = [InPar.Map.Ind]==Iseg;
           if (~isempty(find(MapInd, 1))),

               %SegmentsInfo(Iseg).Map        = InPar.Map(MapInd);
               SegmentsInfo(Iseg).OnSkyX        = InPar.Map(MapInd).X;
               SegmentsInfo(Iseg).OnSkyY        = InPar.Map(MapInd).Y;

               % build the input for griddata/TriScatteredInterp

               All(Iseg).X    = ones(1,Nwave).*SegmentsInfo(Iseg).OnSkyX;
               All(Iseg).Y    = ones(1,Nwave).*SegmentsInfo(Iseg).OnSkyY;
               All(Iseg).Wave = InPar.Wave.';
               All(Iseg).Int  = SegmentsInfo(Iseg).ResampSpec(:,2).';

               %Cube(:,InPar.Map(MapInd).X,InPar.Map(MapInd).Y) = ResampSpec;
           end
       end
    end
end

% build the cube
X    = [All.X].';
Y    = [All.Y].';
Wave = [All.Wave].';
Int  = [All.Int].';


Inn  = ~isnan(X) & ~isnan(Y) & ~isnan(Wave) & ~isnan(Int);
X    = X(Inn);
Y    = Y(Inn);
Wave = Wave(Inn);
Int  = Int(Inn);

[MatX,MatY] = meshgrid(InPar.VecOnSkyX,InPar.VecOnSkyY);

% Back image
Xb = [SegmentsInfo.OnSkyX].';
Yb = [SegmentsInfo.OnSkyY].';
Bb = [SegmentsInfo.MedBack].';
FlagB = ~isnan(Xb) & ~isnan(Yb) & ~isnan(Bb);

Iw = find(Wave==InPar.Wave(100));
switch lower(InPar.CubeInterp)
     case 'griddata'           
          Cube.Back = griddata(Xb(FlagB),...
                               Yb(FlagB),...
                               Bb(FlagB),...
                               MatX,MatY,InPar.Interp2Method);
     case 'tri'
          F = TriScatteredInterp(Xb(FlagB),...
                               Yb(FlagB),...
                               Bb(FlagB));
          Cube.Back = F(MatX,MatY);
     otherwise
            error('Unknown CubeInterp option');
end

Cube.Cube = zeros(Nwave,length(InPar.VecOnSkyY),length(InPar.VecOnSkyX)).*NaN;
Cube.Wave = InPar.Wave;
Cube.Xas  = InPar.VecOnSkyX;
Cube.Yas  = InPar.VecOnSkyY;
% build the cube plane for each wavelength
for Iwave=1:1:Nwave,
    Iw = find(Wave==InPar.Wave(Iwave));
    switch lower(InPar.CubeInterp)
        case 'griddata'           
            Cube.Cube(Iwave,:,:) = griddata(X(Iw),Y(Iw),Int(Iw),MatX,MatY,InPar.Interp2Method);
        case 'tri'
            F = TriScatteredInterp(X(Iw),Y(Iw),Int(Iw));
            Cube.Cube(Iwave,:) = F(MatX,MatY);
        otherwise
            error('Unknown CubeInterp option');
    end
end


% save a FITS cube
if (~isempty(InPar.OutFile)),
    Header = cell_fitshead_addkey([],'COMMENT','','Created by sedm_segments2cube.m written by Eran O. Ofek');
    % add some info from: Cube.KeywordS
    % and WCS
    fitswrite_nd(permute(Cube.Cube,[2 3 1]),InPar.OutFile,Header);
end
