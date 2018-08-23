function [AstC,StFlag]=select_ccdsec(AstC,CCDSEC,Cols)
%--------------------------------------------------------------------------
% select_ccdsec function                                     class/@AstCat
% Description: Select rows from an AstCat object whic are within a
%              CCDSEC region.
% Input  : - An AstCat object.
%          - A CCDSEC vector, or column matrix [Xmin Xmax Ymin Ymax].
%          - Cell array of column names or vector of column indices
%            of the X/Y columns. Default is [1 2].
% Output : - An AstCat object with rows within the CCDSEC coordinates
%            region.
%          - A structure array (elemet per AstCat element) with a field
%            Flag containing the logical flags indicating if the rows
%            in the input AstCat object are within the CCSEC region.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [AstC,StFlag]=select_ccdsec(AstC,[1 100 1 100],{'XWIN_IMAGE','YWIN_IMAGE'});
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<3),
    Cols = [1 2];
end

Nsec = size(CCDSEC,1);

Ncat = numel(AstC);
StFlag = struct('Flag',cell(Ncat,1));
for Icat=1:1:Ncat,
    ColInd = colname2ind(AstC(Icat),Cols);
    Isec   = min(Icat,Nsec);
    
    StFlag(Icat).Flag = AstC(Icat).Cat(:,ColInd(1))>CCDSEC(Isec,1) & ...
                        AstC(Icat).Cat(:,ColInd(1))<CCDSEC(Isec,2) & ...
                        AstC(Icat).Cat(:,ColInd(2))>CCDSEC(Isec,3) & ...
                        AstC(Icat).Cat(:,ColInd(2))<CCDSEC(Isec,4);
    AstC(Icat).Cat = AstC(Icat).Cat(StFlag(Icat).Flag,:);
end

    
    