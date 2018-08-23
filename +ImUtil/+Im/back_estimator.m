function [Mode,Sigma]=back_estimator(Im,varargin)
% Estimate an image background level (mode-like) and noise (std).
% Package: ImUtil.Im
% Description: Given an image, estimate its background level and std by
%              fitting the histogram of the image with a Gaussian. The
%              background is set to the Gaussian peak, while the noise is
%              set to the Gaussian sigma.
% Input  : - A 2-D image.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Buffer' - Number of pixels near image edges that will not be
%                     used in the calculatiom. If empty, set to 0.
%                     Default is empty.
%            'RemUpperQuantile' - Remove upper quantile of pixels before
%                     estimation. Default is 0.9.
%            'MinNumInBin' - Desired mean number of events in bin.
%                     Default is 100.
%                     The final bin size is
%                     max(0.5*ErrOnMean,BinSize_MinNumInBin), where
%                     ErrOnMean is the estimated error on the mean
%                     resolution.
%            'MinRelToMedian' - This parameter multiplied by the image
%                     median, sets the lower range of the pixels histogram.
%                     Default is 0.2.
%            'Nsigma' - This parameter multiplied by the image std plus the
%                     image median, sets the upper range of the pixels
%                     histogram.
%                     Default is 2.
%            'RemoveHighestN' - Number of bins with highest counts to
%                     remove prior to fitting.
%                     Default is 2.
%                     This is useful in order to remove bins with large
%                     anomalous numbers of events due to flat fielding
%                     issues.
%            'RemoveLowN' - Bins with number of events lower than this
%                     threshold will not be used in the fitting.
%                     Default is 10.
%            'Method' - Fitting method:
%                     'polyfit' - Linear fitting of parabola to Gaussian
%                             peak and lines to regions of intersection
%                             with Gaussian sigma.
%                     'fminsearch' - Non linear Gaussian fitting to the
%                             histogram using fminsearch.
%                     Default is 'polyfit'.
%            'OnlyLow' - If method is 'polyfit', this parameter tells if to
%                     use only the low-end part of the Gaussian to estimate
%                     sigma (i.e., mor robust to sources).
%                     Default is false.
% Output : - A robust estimator for the image background level (mode).
%          - A robus estimator for the image noise level (std).
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Mode,Sigma]=ImUtil.Im.back_estimator(Im);
% Reliable: 2
%--------------------------------------------------------------------------


DefV.Buffer               = [];
DefV.RemUpperQuantile     = 0.95;
DefV.MinNumInBin          = 100;
DefV.MinRelToMedian       = 0.2;
DefV.Nsigma               = 2;
DefV.RemoveHighestN       = 2;
DefV.RemoveLowN           = 10;
DefV.Method               = 'fminsearch'; %'polyfit'; %'fminsearch'; %'polyfit';
DefV.OnlyLow              = false;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% remove image edges if needed
if (~isempty(InPar.Buffer))
    Im = Im(InPar.Buffer:end-InPar.Buffer,InPar.Buffer:end-InPar.Buffer);
end

% convert image to vector
Im = Im(:);

if (~isempty(InPar.RemUpperQuantile))
    Th = quantile(Im,InPar.RemUpperQuantile);
    Im = Im(Im<Th);
end

% number of pixels in image
Npix = numel(Im);

% median
Median = median(Im);

%Std    = std(Im);
% faster:
Std   = sqrt(sum((Im - mean(Im)).^2)./(Npix-1));

% best possible knowledge of mode is bounded by the error on the mean
ErrOnMean = Std/sqrt(Npix);

% rough number of binsize to contain MinNumInBin events
%MinBin    = Median./(Npix./InPar.MinNumInBin);
MinBin    = min(2.*Std,Median)./(Npix./InPar.MinNumInBin);

% final bin size
% the 0.5 is for oversampling
BinSize   = max(0.5.*ErrOnMean,MinBin);
%[ErrOnMean,BinSize]

% lower bound on histogram range
HistLowBound = InPar.MinRelToMedian*Median;

% upper bound on histogram range
HistUpperBound = Median + InPar.Nsigma.*Std;

% histogram of counts
%Xbin      = (InPar.RangeRelToMedian(1)*Median:BinSize:InPar.RangeRelToMedian(2)*Median).';
Xbin      = (HistLowBound:BinSize:HistUpperBound).';
Nbin      = histcounts(Im,Xbin).';

Xbin      = (Xbin(1:end-1) + Xbin(2:end)).*0.5;

% remove outliers from histogram
if (InPar.RemoveHighestN>0)
    [~,SI]    = sort(Nbin);
    SI = SI(1:end-InPar.RemoveHighestN);
    Xbin = Xbin(SI);
    Nbin = Nbin(SI);
end

if (InPar.RemoveLowN>0)
    Flag = Nbin>InPar.RemoveLowN;
    Nbin = Nbin(Flag);
    Xbin = Xbin(Flag);
end
    


%[Xbin,Nbin]


switch lower(InPar.Method)
    case 'polyfit'
        % subtract median for better scaling (numerical stability)...
        Xbin = Xbin - Median;
        
        % 2nd order polynomial fit
        % The 0.6 (=exp(-1/2)) is required in order to ensure the 2nd deg polynomial
        % will be around peak and have a maximum (negative 2nd derivative)
        Flag     = Nbin>(0.6.*max(Nbin));
        
        % this is much faster than polyfit
        H = [Xbin(Flag).^2, Xbin(Flag), ones(sum(Flag),1)];
        Par = H\Nbin(Flag);
        %Par      = polyfit(Xbin(Flag),Nbin(Flag),2);
        %Ybin     = H*Par;
        % find polynomial max
        %[~,MaxI] = max(Ybin);
        %Mode     = Xbin(MaxI);
        %MaxVal   = Ybin(MaxI);
        Mode     = Median-0.5.*Par(2)./Par(1);
        MaxVal   = polyval(Par,Mode-Median);
        
        if (nargout>1)
            % estimate sigma
            FlagLow  = Xbin<(Mode-Median);
            ExpHalf  = exp(-0.5);
            FlagSigma = (Nbin./MaxVal)>(ExpHalf-0.2) & (Nbin./MaxVal)<(ExpHalf+0.2);
            FlagLowSigma = FlagSigma & FlagLow;
            Hlow   = [Xbin(FlagLowSigma), ones(sum(FlagLowSigma),1)];
            ParLow = Hlow\Nbin(FlagLowSigma);
            %ParLow = polyfit(Xbin(FlagLowSigma),Nbin(FlagLowSigma),1);
            SigLow = abs((ExpHalf.*MaxVal-ParLow(2))./ParLow(1));

            if (InPar.OnlyLow)
                Sigma = SigLow;
            else
                FlagUp   = Xbin>(Mode-Median);
                if (any(FlagUp))
                    FlagUpSigma  = FlagSigma & FlagUp;
                    Hup    = [Xbin(FlagUpSigma), ones(sum(FlagUpSigma),1)];
                    ParUp  = Hup\Nbin(FlagUpSigma);
                    %ParUp  = polyfit(Xbin(FlagUpSigma),Nbin(FlagUpSigma),1);
                    SigUp  = abs((ExpHalf.*MaxVal-ParUp(2))./ParUp(1)  );
                    Sigma  = 0.5.*(SigLow+SigUp);
                else
                    Sigma = SigLow;
                end
            end
        end
        
%          if (isinf(Sigma))
%              'h'
%          end
       
%         plot(Xbin,Nbin./MaxVal,'.');
%         hold on;
%         plot(Xbin,Ybin./MaxVal,'.');
        

    case 'fminsearch'
        % sort - otherwise will mess up amplitude estimate
        [Xbin,SI] = sort(Xbin);
        Nbin      = Nbin(SI);
        Res=Util.fit.fit_gauss1d(Xbin,Nbin,1);
        if (Res.X0<min(Xbin) || Res.X0>max(Xbin))
            warning('Best fit mode is out of range - set mode to max value and std to sqrt(max)');
            [~,MaxI] = max(Nbin);
            Mode   = Xbin(MaxI);
            Sigma  = sqrt(Mode);
        else
            %plot(Res.Resid)
            %input('hi')
            Mode   = Res.X0;
            Sigma  = Res.Sigma;
        end
    otherwise
        error('Unknown Method option');
end
        
        





