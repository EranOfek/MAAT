function [Sim,Res]=psf_extractor(Sim,varargin)
% Estimate PSF for SIM images.
% Package: @SIM
% Description: Extract the PSF from an astronomical images stored in SIM
%              object. The function first select PSF stars candidates,
%              collect a cube of these stars and construct the mean PSF and
%              its error.
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'StampHalfSize' - PSF Stamp half size. Default is 10
%                              (i.e., PSF size is 21x21).
%            'PosCols' - X/Y Column names by which to extract sources
%                        positions. Default is {'XWIN_IMAGE','YWIN_IMAGE'}.
%            'ExtractorFun' - Source extractor function.
%                        Default is @mextractor.
%            'ExtractorPar' - Cell array of additional arguments to pass
%                        to the source extraction function.
%                        Default is
%                        {'FilterFind',false,'CleanByBackGrad',false,'FlagDiffSpike',false}.
%            'PSFSelectorFun' - Function handle for PSF selection.
%                        Default is @psf_cat_selector.
%            'PSFSelectorPar'- Cell array of additional arguments to pass
%                        to the PSF source selector function.
%                        Default is {}.
%            'PSFCombinerPar' - Cell array of additional arguments to pass
%                        to the PSF combiner function.
%                        Default is {}.
% Output : - A Sim object with the PSF populated.
%          - A structure array with information from the PSF combining
%            process.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=psf_extractor(Sim);
% Reliable: 2
%--------------------------------------------------------------------------


CatField = AstCat.CatField;

DefV.StampHalfSize        = 10;
DefV.PosCols              = {'XWIN_IMAGE','YWIN_IMAGE'};
DefV.ExtractorFun         = @mextractor;
DefV.ExtractorPar         = {'FilterFind',false,'CleanByBackGrad',false,'FlagDiffSpike',false};
DefV.PSFSelectorFun       = @psf_cat_selector;
DefV.PSFSelectorPar       = {};
DefV.PSFCombinerPar       = {};

InPar = InArg.populate_keyval(DefV,varargin,mfilename);



% Find stars (if not exist)
IsCatPop = isfield_populated(Sim,CatField);
if (~all(IsCatPop))
    Sim(~IsCatPop) = InPar.ExtractorFun(Sim(~IsCatPop),InPar.ExtractorPar{:});
end

% for each image
Nsim = numel(Sim);
if (Nsim==0)
    Res = [];
end
for Isim=1:1:Nsim
    % select best stars
    [AstCatP,~] = InPar.PSFSelectorFun(Sim(Isim),InPar.PSFSelectorPar{:});

    % prepare cube of stars
    XY = col_get(AstCatP,InPar.PosCols);
    CubePSF = stamp_xy(Sim(Isim),XY,'StampSize',InPar.StampHalfSize,'Align',false,'HalfPixShift',false,'OutType','cube','UpdateWCS',false);

    % combine cube to PSF
    [CPSF,Res(Isim)]=ImUtil.Im.psf_combiner(CubePSF,InPar.PSFCombinerPar{:});

    % add PSF to Sim
    Sim(Isim).PSF    = CPSF.PSF;
    Sim(Isim).ErrPSF = CPSF.ErrPSF;
end
