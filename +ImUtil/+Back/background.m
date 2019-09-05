function [BackImage,NoiseImage]=background(Image,varargin)
% Calculate the background and noise images for an image
% Package: ImUtil.Back
% Description: Calculate the background and noise images for an image.
% Input  : - A 2D image.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'BackMethod'
%            'NoiseMethod'
% Output : - Background image.
%          - Noise image.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [BackImage,NoiseImage]=ImUtil.Back.background(Image);
% Reliable: 
%--------------------------------------------------------------------------


DefV.BackMethod           = 'median';
DefV.NoiseMethod          = 'std';
DefV.FunOnNoise           = [];  % e.g., @sqrt
DefV.BlockSize            = [128 128];   % [X,Y]   - [] for global
DefV.Buffer               = 10;

DefV.BackQuantilePar      = 0.2;
DefV.BackGaussFitPercent  = [0.025 0.9];
DefV.BackRMeanP           = [0.05 0.05];
DefV.BackNinBin           = 30;

DefV.InterpGenMethod      = 'interp';
DefV.InterpMethod         = 'makima';
DefV.MinNoise             = 1;   % if not empty

InPar = InArg.populate_keyval(DefV,varargin,mfilename);


SizeImage = fliplr(size(Image));  % [X,Y]

% set BlockSize if empty
if isempty(InPar.BlockSize)
    % calculate global value
    InPar.BlockSize = SizeImage;
end


for Icalc=1:1:2
    % two calculations
    % 1 for background
    % 1 for noise
    
    switch lower(InPar.BackMethod)
        case 'mean'
            Fun       = @nanmean;
            FunAddPar = {'all'};
        case 'median'
            Fun       = @nanmedian;
            FunAddPar = {'all'};     
        case 'rmean'
            Fun       = @Util.stat.rmean;
            FunAddPar = {'all',InPar.BackRMeanP};
        case 'std'
            Fun       = @nanstd;
            FunAddPar = {0,'all'};
        case 'var'
            Fun       = @nanvar;
            FunAddPar = {0,'all'};    
        case 'min'
            Fun       = @nanmin;
            FunAddPar = {[],'all'};
        case 'max'
            Fun       = @nanmax;
            FunAddPar = {[],'all'};
        case 'range'
            Fun       = @range;
            FunAddPar = {'all'};
        case 'quantile'
            Fun       = @quantile;
            FunAddPar = {InPar.BackQuantilePar,'all'};
        case 'iqr'
            % inter quantile range 0.25 to 0.75
            Fun       = @iqr;
            FunAddPar = {'all'};
            
        case 'medabsdev'
            % median absolute deviation
            Fun       = @mad;
            FunAddPar = {1,'all'};
        case 'meanabsdev'
            % mean absolute deviation
            Fun       = @mad;
            FunAddPar = {0,'all'};

        case 'fitgaussmode'
            % note that this fun return [Mode, Std]

            Fun        = @Util.stat.mode_fit;
            FunAddPar  = {'TrimEdge2D',0,'Percent',InPar.BackGaussFitPercent,'JoinOut',true};

        case 'mode_density'
            % note that this fun return [Mode, Std]
            
            Fun        = @Util.stat.mode_density;
            FunAddPar  = {'NinBin',InPar.BackNinBin,'JoinOut',true};

        otherwise
            error('Unknown BackMethod option');
    end
    if (Icalc==1)
        BackFun       = Fun;
        BackFunAddPar = FunAddPar;
    else
        NoiseFun       = Fun;
        NoiseFunAddPar = FunAddPar;
    end
end
        
% image partitioning
[Sub,ListEdge,ListCenter] = ImUtil.Back.image_partitioning(Image,InPar.BlockSize,'Buffer',InPar.Buffer);

% apply function
switch lower(InPar.BackMethod)
    case {'fitgaussmode','mode_density'}
        % treat functions in which the std is copled to the mean value
        % (regardless of the user NoiseMethod
        % these are functions that returns two valies [mean, std]
        
        [Val,X,Y,SubFun,Sub]      = ImUtil.Back.image_partitioning_fun(Sub,InPar.BlockSize,BackFun,'FunAddPar',BackFunAddPar,'Buffer',InPar.Buffer);

        
    otherwise
        % treat functions that return one value
        % in this case a seperate call for mean and std is exaceuted
        [Val,X,Y,SubFun,Sub] = ImUtil.Back.image_partitioning_fun(Sub,InPar.BlockSize,BackFun,'FunAddPar',BackFunAddPar,'Buffer',InPar.Buffer);
        [Val(2),~,~,~,~]     = ImUtil.Back.image_partitioning_fun(Sub,InPar.BlockSize,NoiseFun,'FunAddPar',NoiseFunAddPar,'Buffer',InPar.Buffer);
end

        
[BackImage]  = ImUtil.Back.interp_sparse2full(X,Y,Val(1).Val,SizeImage,'Method',InPar.InterpGenMethod,'InterpMethod',InPar.InterpMethod);
[NoiseImage] = ImUtil.Back.interp_sparse2full(X,Y,Val(2).Val,SizeImage,'Method',InPar.InterpGenMethod,'InterpMethod',InPar.InterpMethod);

if ~isempty(InPar.FunOnNoise)
    NoiseImage = InPar.FunOnNoise(NoiseImage);
end

if ~isempty(InPar.MinNoise)
    NoiseImage(NoiseImage<InPar.MinNoise) = InPar.MinNoise;
end

