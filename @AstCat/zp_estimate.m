function [MagZP,FluxZP,Res]=zp_estimate(AstC,IndRef,varargin)
% Estimate the relative photometric zeropoints for catalogs.
% Package: @AstCat
% Description: Estimate the relative photometric zeropoints for catalogs.
% Input  : - AstCat object containing multiple catalogs for which relative
%            zeropoints will be calculated.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'SkipWCS' - SkipWCS option in the AstCat/match.
%                        Default is true.
%            'MatchPar'- Cell array of additional arguments to pass to
%                        AstCat/match.
%                        Default is {}.
% Output : - Reltaive magnitude ZP.
%          - Relative flux ZP.
%            Smaller number refers to lower transperancy.
%          - Structure of additional information.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [MagZP,FluxZP,Res]=zp_estimate(AstC,IndRef);
% Reliable: 2
%--------------------------------------------------------------------------

DefV.SkipWCS              = true;
DefV.MatchPar             = {};
DefV.Columns              = {'XWIN_IMAGE','YWIN_IMAGE','MAG_PSF','MAGERR_PSF'};
DefV.MagCol               = 'MAG_PSF';
DefV.MagErrCol            = 'MAGERR_PSF';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% Match catalogs
[AstOut,~] = match(AstC,AstC(IndRef),'SkipWCS',InPar.SkipWCS,InPar.MatchPar{:});

[MatchedMat,Summary,~] = astcat2matched_array(AstOut,InPar.Columns);

% Estimahe relative zero point
MagZP = zeros(Summary.Nepoch,1);  % ZP [mag[ relative to IndRef image
Res   = Util.struct.struct_def({'StdDiffMag'},Summary.Nepoch,1);
for Iep=1:1:Summary.Nepoch
    Diff    = MatchedMat.(InPar.MagCol)(:,Iep) - MatchedMat.(InPar.MagCol)(:,IndRef);
    DiffErr = sqrt(MatchedMat.(InPar.MagErrCol)(:,Iep).^2 + MatchedMat.(InPar.MagErrCol)(:,IndRef).^2);
    NN      = ~isnan(Diff);
    Diff    = Diff(NN);
    DiffErr = DiffErr(NN);

    %ZP(Iep) = nanmedian(Diff)
    %MagZP(Iep) = Util.stat.wmedian(Diff,DiffErr);
    MagZP(Iep) = median(Diff);
    Res(Iep).StdDiffMag = nanstd(Diff);
end
% relative ZP in flux units
FluxZP = 10.^(-0.4.*MagZP);