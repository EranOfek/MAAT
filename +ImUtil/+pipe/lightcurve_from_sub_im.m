function [DataT,DataR]=lightcurve_from_sub_im(TargetXY,S,D,Scorr,SigmaF,Summary,CoaddSim,GradSN,varargin)
% Generate light curve from subtraction images
% Package: ImUtil
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% Data=ImUtil.pipe.lightcurve_from_sub_im(TargetXY,S,D,Scorr,SigmaF,CoaddSim);
% Reliable: 
%--------------------------------------------------------------------------

ImageField = SIM.ImageField;

%input : S, D, Scorr, SigmaF, Summary, CoaddSim, TargetXY=[TargetX, TargetY]

DefV.ZP                   = 22;
DefV.SigmaX               = 0;
DefV.SigmaY               = 0;
DefV.Nrand_xy             = 1000;
DefV.RandEdgeBuffer       = 20;
DefV.SubBackS             = false;
DefV.SrcWin               = 2;  % window size for maximum value calc
DefV.MomRadius            = 2;  % window size for moment calculation
DefV.MomSigma             = 2;
DefV.MomMaxIter           = 3;
DefV.MomDetThresh         = 5;
DefV.MinMomThreshPhot     = 3;
DefV.NsigDev0             = 2.5;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);





% store the image subtraction products
% ResSim(Iband).D        = D;
% ResSim(Iband).S        = S;
% ResSim(Iband).Scorr    = Scorr;
% ResSim(Iband).SigmaF   = SigmaF;
% ResSim(Iband).Summary  = Summary;
% ResSim(Iband).CoaddSim = CoaddSim;

% subtract residual background from S
if (InPar.SubBackS)
    % calculate S background
    % suppose to be 0
    S = background(S);
    S = sub_background(S);
end
        
% prepare list of positions at which to measure light curves
SizeIm = imagesize(CoaddSim);  % size [X,Y]
Nim = numel(S);

%TargetXY = 
RandXY   = rand(InPar.Nrand_xy,2).*(SizeIm-2.*InPar.RandEdgeBuffer)+InPar.RandEdgeBuffer;
XY       = [TargetXY; RandXY];
Ntarget  = size(TargetXY,1);
%SourceXY = [ones(size(TargetXY,1),1); 2.*ones(InPar.Nrand_xy,1)];  % 1 - for user; 2 - random

%RoundXY     = round(XY);

Data.X        = XY(:,1);
Data.Y        = XY(:,2);
%Data.SourceXY = SourceXY;
%Nsrc          = size(Data.TargetXY,1);


% read the valus of S at the transient location and normalize ...

% measure the source flux via the S image
% use 


F_S     = [Summary.F_S];
F_S     = F_S(:);

% Random positions
DataR.X          = RandXY(:,1);
DataR.Y          = RandXY(:,2);
DataR.ValTar     = nan(Nim,InPar.Nrand_xy);
DataR.ValErrTar  = nan(Nim,InPar.Nrand_xy);
DataR.ValTarSc   = nan(Nim,InPar.Nrand_xy);
DataR.MaskTar    = zeros(Nim,InPar.Nrand_xy,'uint32');
for Id=1:1:Nim
    [Val,~,MaskVal]      = get_value(S(Id),     [RandXY(:,1), RandXY(:,2)]);
    ValSig               = get_value(SigmaF(Id),[RandXY(:,1), RandXY(:,2)]);
    ValSc                = get_value(Scorr(Id), [RandXY(:,1), RandXY(:,2)]);
    DataR.ValTar(Id,:)    = squeeze(Val).';
    DataR.MaskTar(Id,:)   = squeeze(MaskVal).';
    DataR.ValErrTar(Id,:) = squeeze(ValSig).';
    DataR.ValTarSc(Id,:)  = squeeze(ValSc).';
    
end
DataR.ValTar = DataR.ValTar./F_S;
   
% chi2 for random positions
DataR.Nrand     = InPar.Nrand_xy;
DataR.Nepoch    = numel(F_S);
DataR.EpochChi2 = sum((abs(DataR.ValTar)./DataR.ValErrTar).^2,2);
DataR.PosChi2   = sum((abs(DataR.ValTar)./DataR.ValErrTar).^2,1);



% targets position
DataT.MomX       = nan(Nim,Ntarget);
DataT.MomY       = nan(Nim,Ntarget);
DataT.MomX2      = nan(Nim,Ntarget);
DataT.MomY2      = nan(Nim,Ntarget);
DataT.MomXY      = nan(Nim,Ntarget);
DataT.ValTar     = nan(Nim,Ntarget);
DataT.ValErrTar  = nan(Nim,Ntarget);
DataT.ValTarSc   = nan(Nim,Ntarget);
DataT.MaskTar    = zeros(Nim,Ntarget,'uint32');
DataT.ValMom     = nan(Nim,Ntarget);
DataT.MaskMom    = zeros(Nim,Ntarget,'uint32');
DataT.ValErrMom  = nan(Nim,Ntarget);
for Id=1:1:Nim
    % calculate moments and aperture photometry
    [Mom(Id),Mom2(Id),Aper(Id)]=ImUtil.Im.im_moments(D(Id).(ImageField),TargetXY(:,1),TargetXY(:,2),InPar.MomRadius,InPar.MomSigma,InPar.MomMaxIter);

    %                X Y
    DataT.MomX(Id,:)     = [Mom(Id).X(:)].';
    DataT.MomY(Id,:)     = [Mom(Id).Y(:)].';
    DataT.MomX2(Id,:)    = [Mom2(Id).X2(:)].';
    DataT.MomY2(Id,:)    = [Mom2(Id).Y2(:)].';
    DataT.MomXY(Id,:)    = [Mom2(Id).XY(:)].';
    %Data.MaxX(Id,:)     = MaxX(:).';
    %Data.MaxY(Id,:)     = MaxY(:).';

    % Flux and errors at Target X/Y positions:
    [Val,~,MaskVal]      = get_value(S(Id),     [TargetXY(:,1), TargetXY(:,2)]);
    ValSig               = get_value(SigmaF(Id),[TargetXY(:,1), TargetXY(:,2)]);
    ValSc                = get_value(Scorr(Id), [TargetXY(:,1), TargetXY(:,2)]);
    DataT.ValTar(Id,:)    = squeeze(Val).';
    DataT.MaskTar(Id,:)   = squeeze(MaskVal).';
    DataT.ValErrTar(Id,:) = squeeze(ValSig).';
    DataT.ValTarSc(Id,:)  = squeeze(ValSc).';

    for Isrc=1:1:Ntarget
        if (~isnan(Mom(Id).X(Isrc)) && ~isnan(Mom(Id).Y(Isrc)))
            if (Mom(Id).X(Isrc)>1 && Mom(Id).X(Isrc)<SizeIm(1) &&  Mom(Id).Y(Isrc)>1 && Mom(Id).Y(Isrc)<SizeIm(2))
                [DataT.ValMom(Id,Isrc),~,DataT.MaskMom(Id,Isrc)] = get_value(S(Id),     [DataT.MomX(Id,Isrc), DataT.MomY(Id,Isrc)]);
                DataT.ValErrMom(Id,Isrc)                         = get_value(SigmaF(Id),[DataT.MomX(Id,Isrc), DataT.MomY(Id,Isrc)]);
            end
        end
    end
end
DataT.ValTar = DataT.ValTar./F_S;
DataT.ValMom = DataT.ValMom./F_S;


% add data to DataT:
% Extract value of S at the mean moment position
% for all targets (not random) with Sc>3
Tmp             = DataT.ValTarSc;
FlagNonDet      = Tmp< InPar.MomDetThresh;
Tmp(FlagNonDet) = NaN;
DataT.MeanMomX  = nanmedian(DataT.MomX,1);
DataT.MeanMomY  = nanmedian(DataT.MomY,1);

% Flux and errors at Target X/Y positions:
DataT.ValMeanMom     = nan(Nim,Ntarget);
DataT.MaskMeanMom    = zeros(Nim,Ntarget,'uint32');
DataT.ValErrMeanMom  = nan(Nim,Ntarget);
DataT.ValMeanMomSc   = nan(Nim,Ntarget);
for Id=1:1:Nim
    [Val,~,MaskVal]      = get_value(S(Id),     [DataT.MeanMomX(:), DataT.MeanMomY(:)]);
    ValSig               = get_value(SigmaF(Id),[DataT.MeanMomX(:), DataT.MeanMomY(:)]);
    ValSc                = get_value(Scorr(Id), [DataT.MeanMomX(:), DataT.MeanMomY(:)]);
    DataT.ValMeanMom(Id,:)      = squeeze(Val)./F_S(Id);
    DataT.MaskMeanMom(Id,:)     = squeeze(MaskVal);
    DataT.ValErrMeanMom(Id,:)   = squeeze(ValSig);
    DataT.ValMeanMomSc(Id,:)    = squeeze(ValSc);
    
   
end
DataT.ValBest                = DataT.ValMom;
DataT.ValErrBest             = DataT.ValErrMom;
Dist                         = sqrt((DataT.MomX-DataT.MeanMomX).^2 + (DataT.MomY-DataT.MeanMomY).^2);
FlagNonDet                   = DataT.ValTarSc<InPar.MinMomThreshPhot | Dist>1.5 | isnan(DataT.MomX) | isnan(DataT.MomY);
DataT.ValBest(FlagNonDet)    = DataT.ValTar(FlagNonDet);
DataT.ValErrBest(FlagNonDet) = DataT.ValErrTar(FlagNonDet);

DataT.BestX = nan(Nim,Ntarget);
DataT.BestY = nan(Nim,Ntarget);
DataT.X     = nan(1,Ntarget);
DataT.Y     = nan(1,Ntarget);

for Isrc=1:1:Ntarget
    % for all epochs
    DataT.BestX(:,Isrc)               = DataT.MomX(:,Isrc);
    DataT.BestY(:,Isrc)               = DataT.MomY(:,Isrc);
    DataT.X(Isrc)                     = TargetXY(Isrc,1);
    DataT.Y(Isrc)                     = TargetXY(Isrc,2);
    DataT.BestX(FlagNonDet,Isrc)      = DataT.X(Isrc);
    DataT.BestY(FlagNonDet,Isrc)      = DataT.Y(Isrc);
end

% Search for problematic images
% Images in which the mean random targets is different than zero.
ValTar = DataR.ValTar;        
ValTar(DataR.ValTarSc>3) = NaN;
Tmp    = nanmedian(ValTar,2);

DataT.F_S = F_S;
DataT.RandNon0_Nsig = Tmp./(InPar.NsigDev0.*Util.stat.rstd(Tmp));
DataT.RatioMeanF_S  = DataT.F_S./median(DataT.F_S);

DataT.MagBest    = convert.luptitude(DataT.ValBest,10.^(0.4.*InPar.ZP));
DataT.MagErrBest = 1.086.*DataT.ValErrBest./DataT.ValBest;
   
SGPeak = nan(Nim,Ntarget);
SGx = SIM;
SGy = SIM;
for Id=1:1:Nim
    % propagate astrometric error
    % estimate astrometric noise on photometry
    %[SGx,SGy] = gradient(S(Id));
    SGx.(ImageField) = GradSN(Id).X;
    SGy.(ImageField) = GradSN(Id).Y;
    [SGxPeak] = squeeze(get_value(SGx,[DataT.BestX(Id,:)' DataT.BestY(Id,:)']));
    [SGyPeak] = squeeze(get_value(SGy,[DataT.BestX(Id,:)' DataT.BestY(Id,:)']));
    SGPeak(Id,:) = sqrt(SGxPeak.^2 + SGyPeak.^2);
end

% estimate the phot. err. including astrometric noise and chi2/dof
% correction
SigmaAll        = DataT.ValErrBest.*max(1,sqrt(DataR.EpochChi2./DataR.Nepoch));
SigmaPos        = sqrt(InPar.SigmaX(:).^2 + InPar.SigmaY(:).^2);
SigmaPhotAstrom = SGPeak.*SigmaPos./F_S;
DataT.ValErrAllBest  = sqrt(SigmaAll.^2 + SigmaPhotAstrom.^2);
DataT.MagErrAllBest  = 1.086.*DataT.ValErrAllBest./DataT.ValBest;



% 
%     
%     % write results to file:
%     TextFile = sprintf('LC_Band_%s.txt',Band);
%     FID      = fopen(TextFile,'w');
%     fprintf(FID,'%% Generated by: ImUtil.pipe.imsub_lightcurve\n');
%     fprintf(FID,'%% Eran Ofek\n');
%     fprintf(FID,'%% Generation date: %04d-%02d-%02d %02d:%02d%04.1f\n',clock);
%     fprintf(FID,'%% \n');
%     fprintf(FID,'%% Filter: %s\n',Band);
%     fprintf(FID,'%% Number of images in reference: %d\n',Res(Iband).FluxCoaddNim);
%     fprintf(FID,'%% ZP: %8.3f\n',Res(Iband).ZP);
%     fprintf(FID,'%% JD, Flux, FluxErr, FluxStdRand, Mag, OffsetS, Beta, AM, Az, Alt, Xforce, Yforce, X, Y, X2, Y2, XY, FLAGS\n');
%     fprintf(FID,'%14.6f %10.5f %10.5f %10.5f %8.3f %e %e %5.3f %6.3f %6.3f %7.2f %7.2f %7.2f %7.2f %6.3f %6.3f %6.3f %12d\n',Res(Iband).CatLC.Cat.');
%     fclose(FID);
%     