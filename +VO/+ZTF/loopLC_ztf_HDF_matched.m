function loopLC_ztf_HDF_matched(CatHTM,varargin)
% SHORT DESCRIPTION HERE
% Package: VO
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VO.ZTF.loopLC_ztf_HDF_matched(Cat)
% Reliable: 
%--------------------------------------------------------------------------


DefV.Plot                 = 1;
DefV.MinEpoch             = 5;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


ColCell = {'RA','Dec','I1','I2','Nep','ID','FilterID','Field','RcID','MeanMag','StdMag','RStdMag','MaxMag','MinMag','Chi2','MaxPower','FreqMaxPower'};
Col = cell2struct(num2cell(1:1:numel(ColCell)),ColCell,2);


Nsrc = size(CatHTM,1);
for Isrc=1:1:Nsrc
    
    FieldID = CatHTM(Isrc,Col.Field);
    Lines   = CatHTM(Isrc,[Col.I1, Col.I2]);
    
    [CatLC,ColCellLC,CatProp]=VO.ZTF.read_ztf_HDF_matched(FieldID,Lines);
    
    if (size(CatLC,1)>=InPar.MinEpoch)
        cla
        plot.errorxy(CatLC)

        Isrc
        R = input('any key to cont','s');
    end
    
end
