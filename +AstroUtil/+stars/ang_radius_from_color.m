function [AngRad,Temp]=ang_radius_from_color(MagMat,FamilyCell,BandCell,SystemCell)
% Estimate angular radius and color temperature from a set of magnitudes
% Package: +AstroUtil.stars
% Input  : - A matrix of magnitudes in several bands. Source per line.
%          - A cell array of filter family names for the columns in the
%            marix of magnitudes.
%          - A cell array of band names.
%          - A cell array of mag sys. types.
% Output : - ANgular radius ["]
%          - Temperature [K].
% Example: 
%MagMat=[15.886 16.465 15.168];                                                     
%FamilyCell={'GAIA','GAIA','GAIA'};
%BandCell={'Bp','Rp','G'};
%SystemCell={'Vega','Vega','Vega'};
% [AR,T]=AstroUtil.stars.ang_radius_from_color(MagMat,FamilyCell,BandCell,SystemCell)


RAD = 180./pi;

T  = logspace(log10(3000),log10(30000),100)';
NT = numel(T);

Nband = numel(BandCell);
MagInBand = zeros(NT,Nband);
for Iband=1:1:Nband
    MagInBand(:,Iband) = AstroUtil.spec.blackbody_mag_c(T,FamilyCell{Iband},BandCell{Iband},SystemCell{Iband},constant.SunR,1);
end    
BaseAngeRad = constant.SunR./constant.pc .*RAD.*3600;

Nsrc = size(MagMat,1);
AngRad = zeros(Nsrc,1);
Temp   = zeros(Nsrc,1);
for Isrc=1:1:Nsrc
    RMS = std(MagMat(Isrc,:) - MagInBand,[],2);
    [~,Irms] = min(RMS);
    
    MeanDiff = nanmean(MagMat(Isrc,:) - MagInBand(Irms,:));
    
    AngRad(Isrc) = BaseAngeRad.*10.^(-0.2.*MeanDiff);
    Temp(Isrc)   = T(Irms);
end