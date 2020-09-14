function [Sim,SimBack,SimStD,Info]=background(Sim,varargin)
% Estimate background and std of images in a SIM object.
% Package: @SIM
% Description: Given a SIM object array of images, calculate the local
%              or global background and std for each image, and optionally
%              subtract it from the images.
%              For a simplified version that only subtract the background
%              see sub_background.m.
% Input  : - A SIM object array of images.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField' - The SIM fields on which to execute the operator.
%                          Default is {'Im'}.
%            'SubBack'   - Subtract background from image {true|false}.
%                          Default is false.
%            The following parameters are the input of background_image.m:
%            'MethodBack' - Back level estimator method.
%                       Options are:
%                       'none' - return 0.
%                       'mode_fit' - Use mode_fit.m. 
%                       'gauss_fit'- Use ImUtil.Im.back_estimator. Default.
%                       'mean' - mean
%                       'rmean' - robust mean - see rmean.m for details.
%                       'median' - median
%                       'mode' - Use mode.m
%                       'medfilt2' - 2-D median filter
%                                Pass parameters using 'FunBackPar'.
%                       'ordfilt2' - 2-D order filter
%            'MethodStD' - Std level estimator method.
%                       Options are:
%                       'none' - return 0.
%                       'mode_fit' - Use mode_fit.m.
%                       'gauss_fit'- Use ImUtil.Im.back_estimator. Default.
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
%                      'const'|'triangle'|'dist2','sqrtdist'|'none'.
%                      Default is 'sqrtdist'.
% Output : - The output SIM object array.
%            If SubBack is true then the requested fields (in 'ExecField')
%            will be background subtracted.
%            If only one output argument is requested than the background
%            is populated in the 'BackIm' field and the StD is populated
%            in the 'ErrIm' field.
%          - A SIM object containing the background image in the
%            requested fields (in 'ExecField').
%          - A SIM object containing the StD image in the
%            requested fields (in 'ExecField').
% See also: background_image.m, SIM/sub_background.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=background(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

ErrField    = 'ErrIm';
BackField   = 'BackIm';


DefV.ExecField         = {'Im'};
DefV.SubBack           = false;

DefV.MethodBack         = 'mode_fit';   % 'gauss_fit'|'mode_fit'|'mean'|'median'|'medfilt2'|'ordfilt2'|'mode'
DefV.MethodStD          = 'mode_fit';   % 'gauss_fit'|'mode_fit'|'std'|'rstd'|'sqrt'
DefV.Block              = 'full';
DefV.Buffer             = 0; 
DefV.FunBackPar         = {};
DefV.FunStdPar          = {};
DefV.UseMask            = [];  % SIM, Matrix [1 NaN], Matrix of logicals [true-> use]
DefV.SimField           = 'Im';
DefV.SmoothKernel       = 'sqrtdist'; %'triangle';  % 'const'|'triangle'|'dist2','sqrtdist'|'none'

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end
Nf = numel(InPar.ExecField);

if (nargout>1)
    SimBack = SIM(size(Sim));
    if (nargout>2)
        SimStD = SIM(size(Sim));
    end
end

Nsim = numel(Sim);
for Isim=1:1:Nsim
    % for each SIM element
    for If=1:1:Nf
        % for each SIM field
        % calculate back and std
        [Back,StD,Info]=ImUtil.Im.background_image(Sim(Isim).(InPar.ExecField{If}),'MethodBack',  InPar.MethodBack,...
                                                                         'MethodStD',   InPar.MethodStD,...
                                                                         'Block',       InPar.Block,...
                                                                         'Buffer',      InPar.Buffer,...
                                                                         'FunBackPar',  InPar.FunBackPar,...
                                                                         'FunStdPar',   InPar.FunStdPar,...
                                                                         'UseMask',     InPar.UseMask,...
                                                                         'SimField',    InPar.SimField,...
                                                                         'SmoothKernel',InPar.SmoothKernel);
        % subtract background
        if (InPar.SubBack)
            Sim(Isim).(InPar.ExecField{If}) = Sim(Isim).(InPar.ExecField{If}) - Back;
        end
        % output options
        if (nargout>1)
            SimBack(Isim).(InPar.ExecField{If}) = Back;
            if (nargout>2)
                SimStD(Isim).(InPar.ExecField{If}) = StD;
            else
                Sim(Isim).(ErrField) = StD;
            end
        else
            Sim(Isim).(BackField) = Back;
            Sim(Isim).(ErrField)  = StD;
        end
            
    end
end

                                                                        
        
        
        
        
        
        
    
    
    

