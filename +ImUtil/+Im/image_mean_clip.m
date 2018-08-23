function [Mean,StD,NpixUse]=image_mean_clip(Mat,varargin)
% Sigma clip mean for cube of images.
% Package: ImUtil.Im
% Description: Given a cube of images, calculate the mean image using
%              sigma clip mean. Various sigma clipping methods are
%              supported.
% Input  : - Cube of images. By default the first dimension is the image
%            index (see 'ImDim' option).
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'MeanFun' - Function handle for mean calculation
%                        (e.g., @nanmedian). Default is @nanmean.
%            'RejectMethod' - Clipping rejection method:
%                        'std' - standard deviation
%                        'rstd' - robust std.
%                        'minmax' - min/max number of images.
%                        'perc'   - lower/upper percentile of images.
%            'Reject' - Rejection parameters [low high].
%                       e.g., [3 3] for 'std', [0.05 0.85] for 'perc',
%                       or [1 4] for 'minmax'.
%            'MaxIter'- Number of sigma clipping iterations. 0 for no
%                       sigma clipping. Default is 1.
%            'ImDim'  - The dimension of the image index in the cube.
%                       Default is 1.
% Output : - Mean image
%          - StD image.
%          - Image of number of images used per pixel.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Mean,StD,NpixUse]=ImUtil.Im.image_mean_clip(rand(2,2,100));
% Reliable: 2
%--------------------------------------------------------------------------


DefV.MeanFun          = @nanmean;   % Fun(Mat,Dim)
DefV.RejectMethod     = 'std';      % 'minmax''std','rstd'|'perc',
DefV.Reject           = [3 3];      % low/high rejection
DefV.MaxIter          = 1;          % 0 no iterations
DefV.ImDim            = 1;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Dim = InPar.ImDim;

Mean = (InPar.MeanFun(Mat,Dim));
if (InPar.MaxIter>0 || nargout>1)
    StD  = (std(Mat,[],Dim));
end
Size = size(Mat);
Nim  = Size(Dim);
NpixUse = ones(size(Mean)).*Nim;

switch lower(InPar.RejectMethod)
    case 'perc'
        InPar.Reject = [ceil(InPar.Reject(1).*Nim), ceil(Nim - InPar.Reject(2).*Nim)];
        InPar.RejectMethod = 'minmax';
    otherwise
        % do nothing
end
        
Iter     = 0;
ContLoop = Iter<InPar.MaxIter;
while ContLoop
    Iter = Iter + 1;
    switch lower(InPar.RejectMethod)
        case 'std'
            Flag = bsxfun(@gt,Mat,Mean-StD.*abs(InPar.Reject(1))) & bsxfun(@lt,Mat,Mean+StD*abs(InPar.Reject(2)));
            MatN = Mat;
            MatN(~Flag) = NaN;
            NpixUse = sum(~isnan(MatN),Dim);
            Mean = (InPar.MeanFun(MatN,Dim));
            StD  = (nanstd(MatN,[],Dim));
        case 'rstd'
            Flag = bsxfun(@gt,Mat,Mean-StD.*abs(InPar.Reject(1))) & bsxfun(@lt,Mat,Mean+StD*abs(InPar.Reject(2)));
            MatN = Mat;
            MatN(~Flag) = NaN;
            NpixUse = sum(~isnan(MatN),Dim);
            Mean = (InPar.MeanFun(MatN,Dim));
            StD  = (Util.stat.rstd(MatN,Dim));
        case 'minmax'
            % deal also with 'perc'
            SortedMat = sort(Mat,Dim);
            if (Dim==3)
                SortedMat(:,:,1:InPar.Reject(1))         = NaN;
                SortedMat(:,:,end-InPar.Reject(2)+1:end) = NaN;
            elseif (Dim==1)
                SortedMat(1:InPar.Reject(1),:,:)         = NaN;
                SortedMat(end-InPar.Reject(2)+1:end,:,:) = NaN;
            else
                error('With mimmax option Dim must be either 1 or 3');
            end
            NpixUse = sum(~isnan(SortedMat),Dim);
            Mean = (InPar.MeanFun(SortedMat,Dim));
        case 'none'
            NpixUse = sum(~isnan(Mat),Dim);
            Mean = (InPar.MeanFun(Mat,Dim));
        otherwise
            error('Unknown RejectMethod option');
    end
    ContLoop = Iter<InPar.MaxIter;
end


Mean = squeeze(Mean);
StD  = squeeze(StD);
NpixUse = squeeze(NpixUse);
