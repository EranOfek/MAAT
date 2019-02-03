function [AlSim,ResFit,IndRef]=align(Sim,SimRef,varargin)
% Register a set of SIM images
% Class  : class/@SIM
% Description: Find the best transformation between images and register
%              them. If the transformation between the two images is
%              already known then use SIM/transform.m.
% Input  : - A SIM object.
%          - A SIM object containing a single reference image, or a
%            reference image for each element in the first argument SIM
%            object.
%            If empty then select reference from first SIM input.
%            Default is empty.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'TranC' -    A TranClass object that define the transformation
%                         to fit.
%                         Default is:
%                         TranClass({@FunOne, [],@FunX,[],@FunY,[],@FunTiltXp,[],@FunTiltXn,[]}, {@FunOne, [],@FunX,[],@FunY,[],@FunTiltYp,[],@FunTiltYn,[] });
%            'ExtractorReCalc' - Re execute the source extractor.
%                         Default is false.
%            'ExtractorProg' - Source extractor function handle.
%                         Default is @mextractor.
%            'ExtractorPar' - Cell array of additional arguments to pass to
%                         the source extraction function.
%                         Default is {}.
%            'ChooseRefMethod' - Index of image in SIM to use as a
%                         reference image. Alternatively one of the
%                         following methods:
%                         'max' - Choose the image with the largest number
%                                 of sources.
%                         Default is 'max'.
%            'AlignMethod' - One of the following alignment methods:
%                         'wcsCoo' - Match the sources based on the RA/Dec
%                                    in the catalogs and than fit the X/Y
%                                    transformation.
%                         'xyCoo'  - match the sources based on the X/Y in
%                                    the catalogs and than fit the X/Y
%                                    transformation.
%                                    This option assumes that the catalogs
%                                    are roughly aligned.
%                         'ShiftMatch' - Use the pattern_match_shift.m
%                                    function to find the best shift in X/Y
%                                    and than fit the X/Y transformation.
%                         Default is 'wcsCoo'.
%            'MatchPar' - Cell array of additional arguments to pass to the
%                         AstCat/match function.
%            'MatchedKeys' - Cell array of matched keys (catalog column
%                         names). Default is:
%                         {'XWIN_IMAGE','YWIN_IMAGE','MAG_PSF'}.
%            'FitTransformPar' - Cell array of additional arguments to pass
%                         to the TranClass/fit_transform function.
%                         Default is {}.
%            'ExecField' - Cell array of fields in the SIM object for which
%                         to apply the image transformation.
%                         Default is {'Im'}.
%            'TransformPar' - Cell array of additional arguments to pass to
%                         the SIM/transform.m function.
%                         Default is {}.
%            'TransformImage' - Apply image transformation.
%                         Default is true.
%                         Note that if false than the output image is empty
%                         and only the transformation is calculated.
%            'UpdateCat'- Apply transformation to catalog.
%                         Default is true. NOT IMPLEMENTED YET.
%            'ColX'     - Cell array of X coordinate column names to which
%                         to apply the X transformation.
%                         Default is
%                         {'XWIN_IMAGE','X','XPEAK_IMAGE','X_IMAGE'}.
%            'ColY'     - Cell array of Y coordinate column names to which
%                         to apply the Y transformation.
%                         Default is
%                         {'YWIN_IMAGE','Y','YPEAK_IMAGE','Y_IMAGE'}.
%            'Verbose'  - Verbose. Default is false.
% Output : - A SIM object with the aligned (registered) images.
%          - A structure array containing the best fit transformation
%            information per image.
%            See: TranClass/fit_transform for details.
%          - If relevant - Index of image used as reference image.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=FITS.read2sim('test_PTF_*f02_p004162_c08.fits');
%          S=mextractor(S);
%          [AlSim,ResFit,I]=align(S);
% Reliable: 2
%--------------------------------------------------------------------------


CatField = AstCat.CatField;

DefV.XcooName             = 'XWIN_IMAGE';
DefV.YcooName             = 'YWIN_IMAGE';
DefV.PatMatchPar          = {};

DefV.TranC                = TranClass({@FunOne, [],@FunX,[],@FunY,[],@FunTiltXp,[],@FunTiltXn,[]}, {@FunOne, [],@FunX,[],@FunY,[],@FunTiltYp,[],@FunTiltYn,[] });
DefV.ExtractorReCalc      = false;
DefV.ExtractorProg        = @mextractor;
DefV.ExtractorPar         = {};
DefV.ChooseRefMethod      = 'max';  % or index
DefV.AlignMethod          = 'wcsCoo';  % 'wcsCoo' | 'xyCoo' | 'ShiftMatch'
DefV.MatchPar             = {};
DefV.MatchedKeys          = {'XWIN_IMAGE','YWIN_IMAGE','MAG_PSF'};
DefV.FitTransformPar      = {};
DefV.ExecField            = {SIM.ImageField};
DefV.TransformPar         = {};
DefV.TransformImage       = true;
DefV.UpdateCat            = true;
DefV.ColX                 = {'XWIN_IMAGE','X','XPEAK_IMAGE','X_IMAGE'};
DefV.ColY                 = {'YWIN_IMAGE','Y','YPEAK_IMAGE','Y_IMAGE'};
DefV.Verbose              = false;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (nargin<2)
    SimRef = [];
end
IndRef = [];

Nsim    = numel(Sim);
if (InPar.ExtractorReCalc)
    FlagRecCalcSim    = true(Nsim,1);
else
    FlagRecCalcSim    = ~isfield_populated(Sim,   CatField);
end
% call source extratcor
if (any(FlagRecCalcSim))
    Sim(FlagRecCalcSim)       = InPar.ExtractorProg(Sim(FlagRecCalcSim),      InPar.ExtractorPar{:});
end

% set SimRef
if (isempty(SimRef))
    % choose SimRef from Sim
    if (ischar(InPar.ChooseRefMethod))
        switch lower(InPar.ChooseRefMethod)
            case 'max'
                [~,IndRef] = max(sizecat(Sim));
                SimRef     = Sim(IndRef);
            otherwise
                error('Unknown ChooseRefMethod option');
        end
    else
        % assume ChooseRefMethod is an index
        SimRef = Sim(InPar.ChooseRefMethod);
    end
end


% extract sources in Ref image
NsimRef = numel(SimRef);

if ~(Nsim==NsimRef || NsimRef==1)
    error('Number of elements of SimRef must be 1 or like Sim');
end
if (InPar.ExtractorReCalc)
    FlagRecCalcSimRef    = true(Nsim,1);
else
    FlagRecCalcSimRef    = ~isfield_populated(SimRef,   CatField);
end
% call source extratcor
if (any(FlagRecCalcSimRef))
    SimRef(FlagRecCalcSimRef)       = InPar.ExtractorProg(SimRef(FlagRecCalcSimRef),      InPar.ExtractorPar{:});
end




Nf = numel(InPar.ExecField);  % number of SIM fields on which to execute the transformation


% Matching the stars in the images
switch lower(InPar.AlignMethod)
    case 'wcscoo'
        % Match the stars using their RA/Dec coordinates in catalog
        % this mode relys on the existence of the WCS
        % Note that the WCS is used only in order to match the stars
        % and not inorder to define the transformation

        AstM         = match(Sim,SimRef,InPar.MatchPar{:});
        AstM         = [AstCat.sim2astcat(SimRef); AstM(:)];
        MatchedCat   = astcat2matched_array(AstM,InPar.MatchedKeys);
        ResFit       = fit_transform(InPar.TranC,MatchedCat,InPar.FitTransformPar{:},'SelectRefMethod',1);
        ResFit       = ResFit(2:end);
    case 'xycoo'
        % Match the stars using X/Y coordinates and assuming the catalogs
        % are roughly aligned.
        AstM         = match(Sim,SimRef,InPar.MatchPar{:},'SkipWCS',true);
        AstM         = [AstCat.sim2astcat(SimRef); AstM(:)];
        MatchedCat   = astcat2matched_array(AstM,InPar.MatchedKeys);
        ResFit       = fit_transform(InPar.TranC,MatchedCat,InPar.FitTransformPar{:},'SelectRefMethod',1);
        ResFit       = ResFit(2:end);
        
        
    case 'shiftmatch'
        % Match the stars using the pattern_match_shift function
        % that looks for a consistent shift between the two lists
        % of sources
        
        error('ShiftMatch option is not available yet');
        
%         Nsim = numel(Sim);
%         for Isim=1:1:Nsim
%             Iref     = min(Isim,NsimRef);
%             [Res,H2] = pattern_match_shift(Sim(Isim),SimRef(Iref),InPar.PatMatchPar{:});
%             MatchedCat.(InPar.XcooName) = [Res(1).Res(Res(1).IndBest).MatchedRef(:,1), Res(1).Res(Res(1).IndBest).MatchedCat(:,1)];
%             MatchedCat.(InPar.YcooName) = [Res(1).Res(Res(1).IndBest).MatchedRef(:,2), Res(1).Res(Res(1).IndBest).MatchedCat(:,2)];
%             RF                          = fit_transform(InPar.TranC,MatchedCat,InPar.FitTransformPar{:},'SelectRefMethod',1,'MagField',[]);
%             ResFit(Isim)                = RF(2);
%         end

        

    otherwise
        error('Unknown AlignMethod option');
end

if (InPar.TransformImage)
    
    AlSim = SIM(size(Sim));
    for Isim=1:1:Nsim
        if (InPar.Verbose)
            fprintf('Align image %d\n',Isim);
        end
        % for each image
        %Iref = min(NsimRef,Isim);

        % apply the transformation to the image
        for If=1:1:Nf
            if (isfield_populated(Sim(Isim),InPar.ExecField{If}))
                AlSim(Isim) = transform(Sim(Isim),ResFit(Isim).TranC,InPar.TransformPar{:},'ExecField',InPar.ExecField{If});
            end
        end
        
        % apply the transformation to the catalog
        if (InPar.UpdateCat)
            % apply inverse transformation to catalaog
            AlSim(Isim) = apply_tranclass_inv(AlSim(Isim),ResFit(Isim).TranC,'ColX',InPar.ColX,'ColY',InPar.ColY);
        end

    end
else
    AlSim = [];
end
