function [Peaks,Data]=find_src1d(X,Y,varargin)
%--------------------------------------------------------------------------
% find_src1d function                                               ImSpec
% Description: Find sources (e.g., local maxima above the noise)
%              in 1-D data. Find their location, S/N, width, and
%              local background.
% Input  : - X position vector. If empty then use serial index.
%          - Y vector.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'PeakAlgo' - Peak finding algorithm. One of the following:
%                         'mf' - Matched filter.
%            'Subtract' - Background subtraction algorithm.
%                         See subtract_back1d.m for options.
%                         Default is 'medfilt'
%            'SubPar'   - A cell array of additional parameters to pass
%                         to subtract_back1d.m. Default is {100}.
%            'Thresh'   - Detection threshold in units of background rms.
%            'InterpMethod' - Interpolation method. Default is 'linear'.
%            'MinSep'   - Minimum seperation between peaks. Closer peaks
%                         will be merged. Default is 10.
%            'MedianSm' - median filter smoothing size. Default is 5.
%            'Bias'     - Value of constant bias to subtract before search.
%                         Default is 0.            
%            'GuessWidth' - Gaussian fitting guess width. Default is 5.
%            'BackSize' - Size of background region on each size of the
%                         source. Default is 15.
%            'MaxIntersect' - Maximum number of pixels intersection between
%                         background and sources. Default is 3.
%            'NsigSrc'  - Number of sigma for source region definition.
%                         The sigma referes to the background sigma.
%                         Default is 2.
%            'NsigBck'  - Number of sigma for background region definition.
%                         The sigma referes to the background sigma.
%                         Default is 1.
%            'Gain'     - CCD gain for noise estimatation. Default is 1.
%            'RN'       - CCD readout noise [e-] for noise estimatation.
%                         Default is 10.
%            'MinAperRad'- Minimum aperture radius to select even if
%                         optimal aperture is smaller than this value.
%                         Default is 1.
%            'MaxAperRad'- Maximum aperture radius to test for optimal
%                         aperture. Default is 30.
%            'GoodRangeSpatPos' - Select only peaks within this range
%                         Default is [-Inf Inf].
% Output : - Structure array of peaks information. Each element contains
%            the following fields:
%            .X          - Peak X position.
%            .Y          - Peak Y position.
%            .Ybs        - Peak background subtracted Y position.
%            .D2         - 2nd derivative at peak.
%            .Fit        - Fit structure info returned from fit_specline.m
%            .OptimAperRad - Recomended half size of the source aperture.
%            .DX         - 
%            .BackIndL   - Indices of left background region.
%            .BackIndR   - Indices of right background region.
%            .Back       - 2x2 matrix of background region.
%                          First line for left background region boundries.
%                          Second line for right background region boundries.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Peaks,Data]=find_src1d([],Collpase{1});
% Reliable: 2
%--------------------------------------------------------------------------


DefV.PeakAlgo     = 'rms';
DefV.Subtract     = 'medfilt';
DefV.SubPar       = {100};
%DefV.Filter       = {'fun_gauss',[1 5 2],(1:1:10).',0};
DefV.Thresh       = 8;
DefV.InterpMethod = 'linear';
DefV.MinSep       = 15;
DefV.MedianSm     = 5;
DefV.Bias         = 0;  %786;
%DefV.RN           = 10;
DefV.GuessWidth   = 5;
DefV.BackSize     = 15;
DefV.MaxIntersect = 3;   % maximum number of background pixels allowed to intersect with another source
DefV.NsigSrc      = 2;
DefV.NsigBck      = 0.8;
DefV.Gain         = 1;
DefV.RN           = 10;
DefV.MinAperRad   = 1;
DefV.MaxAperRad   = 30;
DefV.GoodRangeSpatPos = [-Inf Inf];

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

N = numel(Y);
if (isempty(X)),
    X = (1:1:N).';
end

% subtract background
Y       = Y - InPar.Bias;
Ybs     = timeseries.subtract_back1d(Y,InPar.Subtract,InPar.SubPar{:});  % background subtracted
BackY   = Y - Ybs;  % local background estimation
StdY    = stdfilt1(Y);
MedStdY = median(StdY);

%--- find peaks ---
switch lower(InPar.PeakAlgo)
    case 'rms'
        
        if (~isempty(InPar.MedianSm)),
           Ybs    = medfilt1(Ybs,3);
        end
        Extram = find_local_extramum(X,Ybs);
        StdYP  = interp1(X,StdY,Extram(:,1),InPar.InterpMethod);
        
        Extram  = Extram(Extram(:,3)<0 & Extram(:,2)>(StdYP.*abs(InPar.Thresh)),:);
        % remove neighboring peaks
        ExtramC  = zeros(size(Extram));
        FlagDone = zeros(size(Extram,1),1);
        Ind      = 0;
        for Ie=1:1:size(Extram,1),
            Ip = find(abs(Extram(Ie,1) - Extram(:,1))<InPar.MinSep & FlagDone==0);
            FlagDone(Ip) = 1;
            if (~isempty(Ip)),
                if (numel(Ip)>1),
                    [~,IIp] = max(Extram(Ip,2));
                    Ind = Ind + 1;
                    ExtramC(Ind,:) = Extram(Ip(IIp),:);

                else
                    Ind = Ind + 1;
                    ExtramC(Ind,:) = Extram(Ip,:);
                end
            end
        end
        ExtramC = ExtramC(1:Ind,:);  % seperated peaks info [X, Y, 2nd deriv]
        
   
    otherwise
        error('Unknown PeakAlgo option');
end

%--- measure peaks properties ---
Npeak = size(ExtramC,1);
Peaks = struct('X',num2cell(zeros(Npeak,1)),...
               'Y',num2cell(zeros(Npeak,1)),...
               'Ybs',num2cell(zeros(Npeak,1)),...
               'D2',num2cell(zeros(Npeak,1)));

FlagSrc = zeros(size(Y));   % Flag indicating how many source ID in position
for Ipeak=1:1:Npeak,
    Peaks(Ipeak).X   = ExtramC(Ipeak,1);
    Peaks(Ipeak).Y   = interp1(X,Y,Peaks(Ipeak).X);
    Peaks(Ipeak).Ybs = interp1(X,Ybs,Peaks(Ipeak).X);
    Peaks(Ipeak).D2  = ExtramC(Ipeak,3);  % 2nd derivative of peak

    % background range
    
    % fit a gaussian
    Peaks(Ipeak).Fit=fit_specline([X,Ybs],@fun_gauss,'Par0',...
                                  [Peaks(Ipeak).Ybs, Peaks(Ipeak).X, InPar.GuessWidth],...
                                  'Back','none','Plot','n');
    A  = Peaks(Ipeak).Fit.Par(1);  % Gaussian norm
    X0 = Peaks(Ipeak).Fit.Par(2);  % Gaussian center
    S0 = Peaks(Ipeak).Fit.Par(3);  % Gaussian sigma
    
    DX = sqrt(2.*S0.^2.*log(A./(InPar.NsigSrc.*MedStdY)));  % distance from Gaussian center at which Gaussian reach Nsig x noise level
    LineFlag = X>(X0-DX) & X<(X0+DX);
    FlagSrc(LineFlag) = FlagSrc(LineFlag) + 1;
    Peaks(Ipeak).DX = DX;
end


% remove peaks which are out of the good range
GoodPeaks = [Peaks.X]>InPar.GoodRangeSpatPos(1) & [Peaks.X]<InPar.GoodRangeSpatPos(2);
Peaks     = Peaks(GoodPeaks);
Npeak     = numel(Peaks);

PossR = (0:1:InPar.MaxAperRad).';   % possible aperture radius to test

%--- select background regions ---
for Ipeak=1:1:Npeak,
    A  = Peaks(Ipeak).Fit.Par(1);  % Gaussian norm
    X0 = Peaks(Ipeak).Fit.Par(2);  % Gaussian center
    S0 = Peaks(Ipeak).Fit.Par(3);  % Gaussian sigma
    
   
    % distance from Gaussian center at which Gaussian reach Nsig x noise level
    DX = sqrt(2.*S0.^2.*log(A./(InPar.NsigBck.*MedStdY)));  
 
    % First guess for backgroubnd
    Peaks(Ipeak).BackIndL = find(X>(X0 - DX - InPar.BackSize) & X<(X0 - DX));
    Peaks(Ipeak).BackIndR = find(X>(X0 + DX) & X<(X0 + DX + InPar.BackSize));
    
    % check if there are other sources in background region
    
%     IntersectXL = intersect(find(Peaks(Ipeak).BackIndL),find(FlagSrc));
%     if (numel(IntersectXL)>InPar.MaxIntersect),
%         % background overlap with another source
%         % attempt to modify background range
%         Peaks(Ipeak).BackIndL = (max(IntersectXL)+1:1:max(Peaks(Ipeak).BackIndL));
%     end
%     IntersectXR = intersect(find(Peaks(Ipeak).BackIndR),find(FlagSrc));
%     if (numel(IntersectXR)>InPar.MaxIntersect),
%         % background overlap with another source
%         % attempt to modify background range
%         Peaks(Ipeak).BackIndR = (max(IntersectXR)+1:1:max(Peaks(Ipeak).BackIndR));
%     end
    Peaks(Ipeak).Back = [X0 - DX - InPar.BackSize, X0 - DX;    X0 + DX,  X0 + DX + InPar.BackSize];

    
    % calculate optimal aperture (maximize S/N)
    MedianBack = median(Y([Peaks(Ipeak).BackIndR;Peaks(Ipeak).BackIndL]));
    Flux  = A.*(normcdf(PossR,0,S0)-0.5).*2;
    SN    = Flux.*InPar.Gain./sqrt(Flux.*InPar.Gain + 2.*PossR.*MedianBack.*InPar.Gain + 2.*PossR.*InPar.RN.^2);
    %plot(PossR,SN)
    [MaxSN,IndMaxSN] = max(SN);
    Peaks(Ipeak).MaxSN = MaxSN;
    Peaks(Ipeak).OptimAperRad = max(PossR(IndMaxSN),InPar.MinAperRad);
    Peaks(Ipeak).MedianBack   = MedianBack;
end

Ybs_src = Ybs;
Ybs_src(FlagSrc==0) = NaN;
Data.X       = X;
Data.Y       = Y;
Data.Ybs     = Ybs;
Data.Ybs_src = Ybs_src;

% plot(X,Ybs,'k')
% hold on
% Ybs_src = Ybs;
% Ybs_src(FlagSrc==0) = NaN;
% plot(X,Ybs_src,'r-')
% 

