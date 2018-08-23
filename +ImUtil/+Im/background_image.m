function [Back,StD,Info]=background_image(Image,varargin)
% Calculate local background of a 2-D matrix.
% Package: ImUtil.Im
% Description: Given a 2-D image (matrix), calculate its local or global
%              background level and standard deviation (noise) using
%              various estimators.
% Input  : - An image (2-D matrix).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'MethodBack' - Back level estimator method.
%                       Options are:
%                       'none'   - return 0.
%                       'mode_fit' - Use mode_fit.m. Default.
%                       'mean' - mean
%                       'rmean' - robust mean - see rmean.m for details.
%                       'median' - median
%                       'mode' - Use mode.m
%                       'medfilt2' - 2-D median filter
%                                Pass parameters using 'FunBackPar'.
%                       'ordfilt2' - 2-D order filter
%            'MethodStD' - Std level estimator method.
%                       Options are:
%                       'mode_fit' - Use mode_fit.m. Default.
%                       'std' - std
%                       'rstd' - robust std using rstd.m
%                       'sqrt' - Sqrt of back level.
%            'Block' - [X,Y] Block size (sub image size) in which to
%                      calculate the back/std (e.g., [512 512]).
%                      If 'full' then use full image. Default is 'full'.
%                      See ImUtil.Im.image_blocks.m for details.
%            'Buffer'- Buffer size in the boundries of each block which
%                      will be added to the block size in each dimension.
%                      I.e., this is an overlap region between the blocks.
%                      Default is 0.
%                      See ImUtil.Im.image_blocks.m for details.
%            'FunBackPar' - A cell array of additional arguments to pass
%                      to the Background level estimator function.
%                      For example, for 'medfilt2' this should be the
%                      median filter block size (e.g., {[30 30]}).
%                      Default is {}.
%            'FunStdPar' - A cell array of additional arguments to pass
%                      to the std level estimator function.
%                      Default is {}.
%            'UseMask' - An image mask indicating which pixels to use
%                      in the back/std estimation.
%                      This can be:
%                      A matrix image containing either 1 (use pixel) or
%                      NaN (ignore pixel).
%                      A matrix of logicals: true (use), false (ignore).
%                      A SIM image in which contains one of the options
%                      mentioned above.
%                      If empty then use all pixels. Default is [].
%            'SimField' - If 'UseMask' is a SIM this indicate which field
%                      in the SIM to use. Default is 'Im'.
%            'SmoothKernel' - Smoothing kernel to use in the smoothing of
%                      the final back/std image.
%                      Options are:
%                      'const'|'triangle'|'dist2','sqrtdist'|'none'
% Output : - Background matrix.
%          - Noise matrix.
%          - Structure containing additional information.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: % Calculate the sky using mode_fit (default) in 256x256 blocks
%          [Back,StD,Info]=background_image(S(2).Im,'Block',[256 256],'SmoothKernel','none');
%          % Calculate the sky using a median filter (slow)
%          [Back,StD,Info]=background_image(S(2).Im,'Block','full','SmoothKernel','none','MethodBack','medfilt2','FunBackPar',{[30 30]});
% Reliable: 2
%--------------------------------------------------------------------------
import Util.stat.*

DefV.MethodBack         = 'mode_fit';   % 'gauss_fit'|'mode_fit'|'mean'|'median'|'medfilt2'|'ordfilt2'|'mode'
DefV.MethodStD          = 'mode_fit';   % 'gauss_fit'|'mode_fit'|'std'|'rstd'|'sqrt'
DefV.Block              = 'full';
DefV.Buffer             = 0; 
DefV.FunBackPar         = {};
DefV.FunStdPar          = {};
DefV.UseMask            = [];  % SIM, Matrix [1 NaN], Matrix of logicals [true-> use]
DefV.SimField           = 'Im';
DefV.SmoothKernel       = 'triangle';  % 'const'|'triangle'|'dist2','sqrtdist'|'none'

InPar = InArg.populate_keyval(DefV,varargin,mfilename);


ImSize = size(Image);   % image size [I,J]
Size   = fliplr(ImSize);  % image size [X,Y]

if (ischar(InPar.Block))
    switch lower(InPar.Block)
        case 'full'
            ListEdge   = [1 Size(1) 1 Size(2)];    % [Xmin Xmax Ymin Ymax]
            ListCenter = 0.5.*Size;
        otherwise
            error('Unknown Block option');
    end
else
    [ListEdge,ListCenter] = ImUtil.Im.image_blocks(Size,InPar.Block,InPar.Buffer);
end


% UseMask is an image that indicate which pixels should be used
% for the background/std estimation
if (~isempty(InPar.UseMask))
    if (SIM.issim(InPar.UseMask))
        InPar.UseMask = InPar.UseMask.(InPar.SimField);
    end
    
    if (islogical(InPar.UseMask))
        UseMask = ones(size(Image));
        UseMask(~InPar.UseMask) = NaN;
    else
        UseMask = InPar.UseMask;
    end
    
    % put NaNs in pixels that should not be used
    % in the back/std estimation
    Image = Image.*UseMask;
end
    


% for each sub image block
Back    = zeros(ImSize);
StD     = zeros(ImSize);
Nlist   = size(ListEdge,1);
SubBack = zeros(Nlist,1);
SubStd  = zeros(Nlist,1);
for Ilist=1:1:Nlist
    SubImage = Image(ListEdge(Ilist,3):ListEdge(Ilist,4),ListEdge(Ilist,1):ListEdge(Ilist,2));

    % Calculate background
    switch lower(InPar.MethodBack)
        case 'none'
            SubBack(Ilist) = 0;
            CalcStd        = false;
        case 'mode_fit'
            % is mode_fit requested also for the std:
            switch lower(InPar.MethodStD)
                case 'mode_fit'
                    [SubBack(Ilist),SubStd(Ilist)] = Util.stat.mode_fit(SubImage,InPar.FunBackPar{:});
                    CalcStd = false;
                otherwise
                    SubBack(Ilist) = Util.stat.mode_fit(SubImage,InPar.FunBackPar{:});
                    CalcStd = true;
            end
        case 'mode_fit1'
            % is mode_fit requested also for the std:
            switch lower(InPar.MethodStD)
                case 'mode_fit1'
                    %[BackIm,NoiseIm,Back] = ImUtil.Im.background_fit(Image,varargin)
                    %[SubBack(Ilist),SubStd(Ilist)] = Util.stat.mode_fit(SubImage,InPar.FunBackPar{:});
                    CalcStd = false;
                otherwise
                    %SubBack(Ilist) = Util.stat.mode_fit(SubImage,InPar.FunBackPar{:});
                    CalcStd = true;
            end
        case 'gauss_fit'
            switch lower(InPar.MethodStD)
                case 'gauss_fit'
                    [SubBack(Ilist),SubStd(Ilist)] = ImUtil.Im.back_estimator(SubImage,InPar.FunBackPar{:});
                    CalcStd = false;
                otherwise
                    SubBack(Ilist) = ImUtil.Im.back_estimator(SubImage,InPar.FunBackPar{:});
                    CalcStd = true;
            end
        case 'mean'
            SubBack(Ilist) = nanmean(SubImage(:));
            CalcStd = true;
        case 'rmean'
            SubBack(Ilist) = rmean(SubImage(:),1,InPar.FunBackPar{:});
            CalcStd = true;
        case 'median'
            SubBack(Ilist) = nanmedian(SubImage(:));
            CalcStd = true;
        case 'mode'
            SubBack(Ilist) = mode(SubImage(:));
            CalcStd = true;
        case 'medfilt2'
            % This is working only for Block=='full'
            if (Nlist>1)
                error('MethodBack of medfilt2 is working only for a single block');
            end
            SubBack(Ilist) = medfilt2(SubImage,InPar.FunBackPar{:});
            CalcStd = true;
        case 'ordfilt2'
            % This is working only for Block=='full'
            if (Nlist>1)
                error('MethodBack of ordfilt2 is working only for a single block');
            end
            SubBack(Ilist) = ordfilt2(SubImage,InPar.FunBackPar{:});  
            CalcStd = true;
        otherwise
            error('Unknown MethodBack option');
    end
    if (Nlist>1)
        Back(ListEdge(Ilist,3):ListEdge(Ilist,4),ListEdge(Ilist,1):ListEdge(Ilist,2)) = SubBack(Ilist);
    end
    
    if (CalcStd && nargout>1)
        % calculate StD
        switch lower(InPar.MethodStD)
            case 'none'
                SubStd(Ilist) = 0;
            case 'mode_fit'
                [~,SubStd(Ilist)] = Util.stat.mode_fit(SubImage,InPar.FunBackPar{:});
            case 'gauss_fit'
                [~,SubStd(Ilist)] = ImUtil.Im.back_estimator(SubImage,InPar.FunBackPar{:});
            case 'std'
                SubStd(Ilist) = nanstd(SubImage(:));
            case 'rstd'
                SubStd(Ilist) = rstd(SubImage(:));
            case 'sqrt'
                SubStd(Ilist) = sqrt(SubBack(Ilist));
            otherwise
                error('Unknown MethodStD option');
        end
        
    end
    if (nargout>1 && Nlist>1)
        StD(ListEdge(Ilist,3):ListEdge(Ilist,4),ListEdge(Ilist,1):ListEdge(Ilist,2)) = SubStd(Ilist);
    end
end


% merge SubBack and SubBack into an image or scalar.
if (Nlist==1)
    % If Nlist == 1 then return a scalar values for the entire images
    Back = SubBack;
    StD  = SubStd;
else
    if (InPar.Block(1)<0)
        InPar.Block(1) = ImSize(2)./abs(InPar.Block(1));
    end
    if (InPar.Block(2)<0)
        InPar.Block(2) = ImSize(1)./abs(InPar.Block(2));
    end
    VecX = (1:1:InPar.Block(1));
    VecY = (1:1:InPar.Block(2));
    
    [MatX,MatY] = meshgrid(VecX,VecY);
    
    switch lower(InPar.SmoothKernel)
        case {'none','no'}
            Kernel = [];
        case {'triangle','tri'}
            Dist    = sqrt((MatX-0.5.*InPar.Block(1)).^2 + (MatY-0.5.*InPar.Block(2)).^2);
            MaxDist = sqrt((0.5.*InPar.Block(1)).^2 + (0.5.*InPar.Block(2)).^2);
            Kernel  = MaxDist-Dist;
        case 'const'
            Kernel  = ones(size(MatX));
        case 'dist2'
            Dist    = sqrt((MatX-0.5.*InPar.Block(1)).^2 + (MatY-0.5.*InPar.Block(2)).^2);
            MaxDist = sqrt((0.5.*InPar.Block(1)).^2 + (0.5.*InPar.Block(2)).^2);
            Kernel  = (MaxDist-Dist).^2;
        case 'sqrtdist'
            Dist    = sqrt((MatX-0.5.*InPar.Block(1)).^2 + (MatY-0.5.*InPar.Block(2)).^2);
            MaxDist = sqrt((0.5.*InPar.Block(1)).^2 + (0.5.*InPar.Block(2)).^2);
            Kernel  = sqrt(MaxDist-Dist);
        otherwise
            error('Unknown SmoothKernel option');
    end
    
    % If kernel is empty then back/std map are box pixelated
    if (~isempty(Kernel))
        Kernel  = Kernel./sum(Kernel(:));
        Back = Util.external.conv_fft2(Back,Kernel,'reflect');
        StD  = Util.external.conv_fft2(StD,Kernel,'reflect');
    end
    
    
end
 

Info.ListEdge   = ListEdge;
Info.ListCenter = ListCenter;
Info.SubBack    = SubBack;
Info.SubStd     = SubStd;

    


