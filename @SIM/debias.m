function [Sim,BiasSim,BiasOverScan,NotIsBias,Summary]=debias(Sim,BiasSim,varargin)
% Subtract bias image and bias overscan from images.
% Package: @SIM
% Description: Subtract bias image and bias overscan from images.
%              Optionally also create the bias image.
% Input  : - A SIM object.
%          - An optional SIM object containing the bias image.
%            If empty then bias image will be calculated using SIM/bias.m.
%            Default is empty.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'BiasPar'    - Cell array of additional input arguments to
%                           pass to SIM/bias.m. Default is {}.
%            'IsBiasPar'  - Cell array of additional input arguments to
%                           pass to SIM/isbias.m. Default is {}.
%            'NotIsBias'  - A vector of logicals (or indices) indicating
%                           which input SIM images should be bias
%                           subtracted. If empty, then will attempt to
%                           identify the non-bias images and subtract the
%                           bias from these images. If true, then subtract
%                           the bias from all the images.
%                           Default is true.
%            'SubOverScan'-Subtract bias oversan. Default is true.
%            'OverSecPar'- Cell array of additional input arguments to
%                          pass to SIM/bias_overscan.m. Default is {}.
%            'FinalSec'  - The final image section region. This is applied
%                          in bias_overscan.m.
%                          This should be used if the bias overscan region
%                          should be trimmed.
%                          This is either a string containing an header
%                          keyword name containing the bias overscan
%                          region, or [Xmin Xmax Ymin Ymax], or [].
%                          If empty then do not trim image.
%                          Default is 'CCDSEC'.
%            'ReturnAll' - The output Sim contains all the input images
%                          (true) or only the bias subtracted images
%                          (false). Default is true.
% Output : - SIM object with the bias subtracted images.
%          - SIM object with the super bias image.
%          - SIM object with the overscan bias image per each image.
%          - A vector of logocals indicating not bias images (i.e., images
%            from which bias was subtracted).
%          - Summary meta data from bias.m. See bias.m for details.
% See also: bias.m, isbias.m, bias_overscan.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=debias(S,bias(S));
%          S=debias(S);
%          S=debias(S2,[],'FinalSec',[]);
%          % debias P200 LFC images
%          S=FITS.read2sim('ccd*.0.fits');
%          S=trim_image(S,[1 1056 1 2063]); % fix bug in HEADER info
%          [S,Bias,BiasOS,NIB,Sum]=debias(S,[],'IsBiasPar',{'RN',[]},'FinalSec',[1 1024 1 2063]);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1)
    BiasSim = [];
end

DefV.BiasPar            = {};
DefV.IsBiasPar          = {};
DefV.NotIsBias          = true;    % list of images from which to subtract the bias
DefV.SubOverScan        = true;
DefV.OverScanPar        = {};
DefV.FinalSec           = 'CCDSEC';
DefV.ReturnAll          = true;

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% construct super bias
IsBias  = [];
if (isempty(BiasSim))
    IsBias = isbias(Sim,InPar.IsBiasPar{:});
    [BiasSim,Summary] = bias(Sim,'IsBias',IsBias,InPar.BiasPar{:});
else
    Summary = [];
end

% identify images from which to subtract the superbias
if (isempty(InPar.NotIsBias))
    if (isempty(IsBias))
        % IsBias was not populated
        IsBias = isbias(Sim,InPar.IsBiasPar{:});
    end
    NotIsBias = ~IsBias;
elseif (numel(InPar.NotIsBias)==1 && all(InPar.NotIsBias))
    % subtract bias from all the images
    NotIsBias = true(size(Sim));
else
    % Assume NotIsBias is in the correct form
    % containing: vector of logicals or indices
    NotIsBias = InPar.NotIsBias;
end

% subtract super bias
% and combine the mask of the image with that of the bias
% bfun2sim(double(... -> bfun2sim [no need for double conversion]
Sim(NotIsBias) = bfun2sim(Sim(NotIsBias),BiasSim,@minus,'MaskFun',@mask_add,'MaskFunPar',{@bitor});

% subtract bias overscan
BiasOverScan = [];
if (InPar.SubOverScan)
    [Sim(NotIsBias),BiasOverScan] = bias_overscan(Sim(NotIsBias),'FinalSec',InPar.FinalSec,InPar.OverScanPar{:});
end

% return all images or only bias subtracted images
if (~InPar.ReturnAll)
    Sim = Sim(NotIsBias);
end

