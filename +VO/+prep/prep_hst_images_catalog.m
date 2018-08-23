function C=prep_hst_images_catalog(varargin)
% SHORT DESCRIPTION HERE
% Package: VO
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: C=VO.prep.prep_hst_images_catalog
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

CatField = AstCat.CatField;

DefV.DecVec               = (-90:10:90)';
DefV.DirCats              = '/home/eran/matlab/data/+cats/+sources';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Ndec = numel(InPar.DecVec);
for Idec=1:1:(Ndec-1)
    QS = sprintf('SELECT ImageName, RA, Dec, ImageID, Sexcat_NObj FROM Images WHERE Dec>=%f and Dec<%f',InPar.DecVec(Idec),InPar.DecVec(Idec+1));
    
    Tmp = VO.MAST.query_casjobs(QS,'Table','HSCv2','SaveInTable',true,'FormatString','%s %f %f %f %f\n');
    
    if (Idec==1)
        C = Tmp;
    else
        C.(CatField) = [C.(CatField); Tmp.(CatField)];
    end
end

Nim = numel(C.(CatField).ImageName);
SC = regexprep(C.(CatField).ImageName,'"','');
SC = regexp(SC,'_','split');
Camera = cell(Nim,1);
Channel = cell(Nim,1);
for Iim=1:1:Nim
    Camera{Iim} = SC{Iim}{4};
    Channel{Iim} = SC{Iim}{5};
end

UnCam = unique(Camera);
UnCha = unique(Channel);
CamVec = zeros(Nim,1);
ChaVec = zeros(Nim,1);
for I1=1:1:numel(UnCam)
    F = (strcmp(Camera,UnCam{I1}));
    CamVec(F) = I1;
end
for I1=1:1:numel(UnCha)
    F = (strcmp(Channel,UnCha{I1}));
    ChaVec(F) = I1;
end    
 
C.(CatField) = table2array(C.(CatField)(:,2:end));
C.(CatField) = [C.(CatField), CamVec, ChaVec];
C.ColCell    = {'RA','Dec','ImageID','Nsrc','Camera','Channel'};
C = colcell2col(C);
C.(CatField)(:,1:2) = C.(CatField)(:,1:2)./RAD;
C = sortrows(C,2);
C.ColUnits = {'rad','rad','','',''};
C.Name = 'HST images catalog';
C.Source = 'MAST/HSCv2/Images table';
C.UserData.CameraDic  = UnCam;
C.UserData.ChannelDic = UnCha;

% save in +cats dir
cd(InPar.DirCats)
HSTimages = C;

save HSTimages.mat HSTimages

