function [FRA,FDec]=field2coo(VecField,VecCCDID)
% Find the center equatorial J2000 coordinates for PTF fields/CCDIDs
% Package: VO.PTF
% Description: Look for the PTF fields/CCDs coordinates by its ID.
% Input  : - Vector of Field IDs.
%          - Vector of CCD IDs. If empty then will return the Field center.
% Output : - Vector of J2000.0 RA for fields.
%          - Vector of J2000.0 Dec for fields.
% Tested : Matlab 7.13
%     By : Eran O. Ofek                    Sep 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [RA,Dec]=VO.PTF.field2coo(1000,1);
% Reliable: 2
%------------------------------------------------------------------------------

if (nargin==1)
    VecCCDID = [];
end

CatField = AstCat.CatField;
ColField = AstCat.ColField;

% Load the PTF Fields/CCD SDSS coverge into an AstCat object
Fields = AstCat.loadh2astcat('PTF_FieldCCD.hdf5');

Nf = numel(VecField);
FRA  = zeros(size(VecField));
FDec = zeros(size(VecField));
for If=1:1:Nf
    if (isempty(VecCCDID))
        Ind = find(VecField(If)==Fields.(CatField)(:,Fields.(ColField).PTFFIELD));
        RA  = Fields.(CatField)(Ind,Fields.(ColField).RA);
        Dec = Fields.(CatField)(Ind,Fields.(ColField).Dec);
        [CD1,CD2,CD3] = celestial.coo.coo2cosined(RA,Dec);
        [FRA(If),FDec(If)] = celestial.coo.cosined2coo(mean(CD1),mean(CD2),mean(CD3));
    else
        Ind = find(VecField(If)==Fields.(CatField)(:,Fields.(ColField).PTFFIELD) & ...
                   VecCCDID(If)==Fields.(CatField)(:,Fields.(ColField).CCDID));
        FRA(If)  = Fields.(CatField)(Ind,Fields.(ColField).RA);
        FDec(If) = Fields.(CatField)(Ind,Fields.(ColField).Dec);
    end
end

