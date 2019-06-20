function [Res,TranZP,Sim]=phot_zp(Sim,varargin)
% Estimate photometric zero points of SIM images
% Package: class/@SIM
% Description: Estimate photometric zero points of SIM images using
%              relative photometry (Ofek et al. 2011).
% Input  : - SIM images.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExtractorReCalc' - Re-extract sources true|false.
%                           Default is false.
%            'ExtractorProg' - Source extractor program.
%                           Default is @mextractor.
%            'ExtractorPar'  - Cell array of additional parameters to pass
%                           to the source extractor program.
%                           Default is {}.
%            'MatchPar'      - Cell array of additional arguments to pass
%                           to the AstCat/match program.
%            'MatchedFields' - Cell array of fields to generate a match
%                           matrix.
%                           Default is
%                           {'XWIN_IMAGE','YWIN_IMAGE','MAG_PSF','MAGERR_PSF','SN','PEAKF_VALTOT'}.
%            'MagColName' - Name of magnitude column. Default is 'MAG_PSF'.
%            'MagErrColName' - Name of magnitude error column.
%                           Default is 'MAGERR_PSF'.
%            'RelZPPar'   - Additional arguments to pass to the rel_zp
%                           program. Default is {}.
% Output : - The rel_zp.m output.
%          - Relative transperancy of images (i.e., 10.^(-0.4.*Res.ZP)).
%          - The SIM object with populated catalogs.
% Reference: Ofek et al. 2011 (Appendix A).
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

CatField = AstCat.CatField;

DefV.ExtractorReCalc      = false;
DefV.ExtractorProg        = @mextractor;
DefV.ExtractorPar         = {};
DefV.MatchPar             = {};
DefV.MatchedFields        = {'XWIN_IMAGE','YWIN_IMAGE','MAG_PSF','MAGERR_PSF','SN','PEAKF_VALTOT'}; %,'FLUX_PSF'}; %,'FLUX_APER_1','FLUX_APER_2'};
DefV.MagColName           = 'MAG_PSF';
DefV.MagErrColName        = 'MAGERR_PSF';
DefV.RelZPPar             = {};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Nsim = numel(Sim);

% check if catalog field is populated
if (InPar.ExtractorReCalc)
    PopCat = false(Nsim,1);
else
    PopCat = ~isfield_populated(Sim,CatField);
end

% if Cat field is not populated then call mextractor
Sim(PopCat) = InPar.ExtractorProg(Sim(PopCat),InPar.ExtractorPar{:});

% JD of images
JD = julday(Sim);

% match sources by position
[AstOut,AstUM] = match(Sim,[],InPar.MatchPar{:});

InPar.FracAppear  = 1;
InPar.MaxPeakFlux = 20000;

[MatchedArrary,Summary]=astcat2matched_array(AstOut,InPar.MatchedFields);
Flag = sum(MatchedArrary.SN>20,2)>=Summary.Nepoch.*InPar.FracAppear;

%MatchedArrary.PEAKF_VALTOT>InPar.MaxPeakFlux;
MatchedArrary.(InPar.MagColName)(MatchedArrary.PEAKF_VALTOT>InPar.MaxPeakFlux) = NaN;

Res = ImUtil.Im.rel_zp(MatchedArrary.(InPar.MagColName)(Flag,:),MatchedArrary.(InPar.MagErrColName)(Flag,:),InPar.RelZPPar{:});
TranZP = 10.^(-0.4.*Res.ZP);




