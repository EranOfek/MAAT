function [Res,H2]=pattern_match_shift(Cat,Ref,varargin)
% Match two AstCat objects by X/Y shift pattern recognition 
% Package: @AstcCat
% Description: Version of pattern_match_shift.m that works on AstCat
%              objects input.
%              Given two lists containing two dimensional planar
%              coordinates (X, Y), find the possible shift transformations
%              between the two lists, using a histogram of all possible
%              combination of X and Y distances.        
% Input  : - An AstCat object containing a single catalog.
%            Sorted by Y coordinate.
%          - A reference AstCat object containing a single catalog.
%            Sorted by Y coordinate.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExtractorProg' - Source extractor function handle to run in
%                         case the catalog is not populated and the input
%                         is a SIM object.
%                         Default is @mextractor;
%            'ExtractorProgPar' - Cell array of additional arguments to
%                         pass to the source extraction function.
%                         Default is {}.
%            'PatMatchPar' - Cell array of additional arguments to
%                         pass to the ImUtil.Im.pattern_match_shift function.
%                         Default is {}.
%            'ColCatXY' - Cell array of column names or vector of column
%                         indices containing the X/Y columns in the input
%                         catalog. Default is [1 2].
%            'ColRefXY' - Cell array of column names or vector of column
%                         indices containing the X/Y columns in the input
%                         reference. Default is [1 2].
%            'CatSelect'- A string containing a selection criteria for rows
%                         in the input catalog to use in the matching.
%                         E.g., 'APER_MAG>14 & APER_MAG<18'.
%                         Alternatively, this can be a vector of indices or
%                         a vector of logicals indicating which rows to
%                         use. If empty use all. Default is empty.
%            'RefSelect'- A string containing a selection criteria for rows
%                         in the input reference to use in the matching.
%                         E.g., 'APER_MAG>14 & APER_MAG<18'.
%                         Alternatively, this can be a vector of indices or
%                         a vector of logicals indicating which rows to
%                         use. If empty use all. Default is empty.
% Output : - Structure array (element per image) of structure arrays of
%            all the possible shift solutions.
%            The following fields are available:
%          - Structure array (element per image) of the 2D histogram of
%            all possible combination of X and Y distances.         
% See also: ImUtil.Im.pattern_match_shift.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res,IndBest]=pattern_match_shift(Cat,Ref);
% Reliable: 2
%--------------------------------------------------------------------------

CatField         = 'Cat';

DefV.ExtractorProg      = @mextractor;
DefV.ExtractorProgPar   = {};
DefV.PatMatchPar        = {};
DefV.ColCatXY           = [1 2];
DefV.ColRefXY           = [1 2];
DefV.CatSelect          = [];
DefV.RefSelect          = [];
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Ncat = numel(Cat);
Nref = numel(Ref);
if (Nref>1)
    error('Ref should contain a single element');
end



% populate catalogs if needed
if (SIM.issim(Cat))
    % check if Cat catalog is populated
    IsCatCatPopulated = isfield_populated(Cat,CatField);
    % source extraction
    Cat(~IsCatCatPopulated) = InPar.ExtractorProg(Cat(~IsCatCatPopulated),InPar.ExtractorProgPar{:});
end
if (SIM.issim(Ref))
    % check if Ref catalog is populated
    IsRefCatPopulated = isfield_populated(Ref,CatField);
    % source extraction
    Ref(~IsRefCatPopulated) = InPar.ExtractorProg(Ref(~IsRefCatPopulated),InPar.ExtractorProgPar{:});
end
    
% select rows in reference catalog
RefColInd = colname2ind(Ref,InPar.ColRefXY);
if (~isempty(InPar.RefSelect))
    RefFlag   = col_arith(Ref,InPar.RefSelect,'mat');
    Ref       = row_select(Ref,RefFlag);
end

for Icat=1:1:Ncat
    % for each Cat element

    % select rows in input catalog
    CatColInd = colname2ind(Cat(Icat),InPar.ColCatXY);
    if (~isempty(InPar.CatSelect))
        CatFlag   = col_arith(Cat(Icat),InPar.CatSelect,'mat');
        Cat       = row_select(Cat(Icat),CatFlag);
    end

    [Res1,IndBest1,H21] = ImUtil.Im.pattern_match_shift(Cat(Icat).(CatField)(:,CatColInd),...
                                           Ref.(CatField)(:,RefColInd),...
                                           InPar.PatMatchPar{:});
                                       
    Res(Icat).Res = Res1;
    Res(Icat).IndBest = IndBest1;
    if (nargout>1)
        H2(Icat).H2 = H21;
    end
                                       
end










