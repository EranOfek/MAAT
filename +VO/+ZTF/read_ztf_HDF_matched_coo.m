function [Cat,ColCell]=read_ztf_HDF_matched_coo(RA,Dec,Filter,varargin)
% Get ZTF/DR1 light curves for source by coordinates
% Package: VO.ZTF
% Description: Get ZTF/DR1 light curves for source by coordinates
% Input  : - J2000.0 R.A. [radians].
%          - J2000.0 Dec. [radians].
%          - Filter index. Default is 1 (g).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ColCell' - Output columns (cell array).
%            'SearchRadius' - Search radius for source (arcsec).
%                        Default is 2.
% Output : - Light curve matrix.
%          - Cell array of column names
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Cat,ColCell]=VO.ZTF.read_ztf_HDF_matched_coo(1.4015719,0.5171304,1)
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin<3)
    Filter = 1;
end

DefV.ColCell            = {'HMJD','Mag','MagErr','ColorCoef','Flags','Field','RcID','Isrc'};
DefV.SearchRadius       = 2;  % arcsec
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (isempty(Dec))
    Dec = RA(:,2);
    RA  = RA(:,1);
end

[CatSrc,ColCell] = catsHTM.cone_search('ztfSrcLCDR1',RA,Dec,InPar.SearchRadius);
Col = cell2struct(num2cell(1:1:numel(ColCell)),ColCell,2);

% select observations at the requested filter
FlagFilter = CatSrc(:,Col.FilterID)==Filter;
CatSrc     = CatSrc(FlagFilter,:);

Nsrc = size(CatSrc,1);
Ncol = numel(InPar.ColCell);
Cat  = zeros(0,Ncol);
for Isrc=1:1:Nsrc
    CatLC = VO.ZTF.read_ztf_HDF_matched(CatSrc(Isrc,Col.Field),CatSrc(Isrc,[Col.I1, Col.I2]));
    Ns  = size(CatLC,1);
    CatLC = [CatLC, ones(Ns,1).*CatSrc(Isrc,[Col.Field, Col.RcID]), ones(Ns,1).*Isrc];
    Cat   = [Cat; CatLC];
end
ColCell = InPar.ColCell;


