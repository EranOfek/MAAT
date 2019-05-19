function [Out]=wget_ps1_api(RA,Dec,varargin)
% Query the PS1 catalog via the web API
% Package: VO
% Description: Query the PS1 catalog via the web API
%              See API details in: https://catalogs.mast.stsci.edu/docs/panstarrs.html
%              and https://archive.stsci.edu/panstarrs/
% Input  : - J2000.0 R.A. [deg].
%          - J2000.0 Dec. [deg].
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Radius' - Search radius [deg]. Default is 1./60.
%            'ColCell'- Cell array of column names to retieve.
%            'BaseURL'- Query base URL.
%            'DR'     - Data Release. Default is 'dr2'.
%            'Catalog'- PS1 catalog. Default is 'maen'.
%            'OutType'- Output type. Default is 'astcat'.
% Output : - PS1 catalog.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Out= VO.MAST.wget_ps1_api(1,1);
% Reliable: 2
%--------------------------------------------------------------------------



DefV.Radius               = 1./60;  % deg
DefV.ColCell              = {'raMean','decMean','raMeanErr','decMeanErr','epochMean','nDetections','qualityFlag',...
                             'gMeanPSFMag','gMeanPSFStd','gMeanPSFNpt','gMeanPSFMagMin','gMeanPSFMagMax','gMeanKronMag','gMeanKronMagErr','gMeanApMag','gMeanApMagErr',...
                             'rMeanPSFMag','rMeanPSFStd','rMeanPSFNpt','rMeanPSFMagMin','rMeanPSFMagMax','rMeanKronMag','rMeanKronMagErr','rMeanApMag','rMeanApMagErr',...
                             'iMeanPSFMag','iMeanPSFStd','iMeanPSFNpt','iMeanPSFMagMin','iMeanPSFMagMax','iMeanKronMag','iMeanKronMagErr','iMeanApMag','iMeanApMagErr',...
                             'zMeanPSFMag','zMeanPSFStd','zMeanPSFNpt','zMeanPSFMagMin','zMeanPSFMagMax','zMeanKronMag','zMeanKronMagErr','zMeanApMag','zMeanApMagErr',...
                             'yMeanPSFMag','yMeanPSFStd','yMeanPSFNpt','yMeanPSFMagMin','yMeanPSFMagMax','yMeanKronMag','yMeanKronMagErr','yMeanApMag','yMeanApMagErr'};
DefV.BaseURL              = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/';
DefV.DR                   = 'dr2';
DefV.Catalog              = 'mean';  % 'mean' | 'stack' | ...
DefV.nDetections          = 2; % 2 or more detections
DefV.pagesize             = '500001';
DefV.TimeOut              = 60; %s
DefV.Format               = 'json'; %'csv';
DefV.OutType              = 'astcat'; % 'astcat' | 'json' | 'mat'
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

ColCell = InPar.ColCell;
Ncol    = numel(ColCell);
Columns = ColCell{1};
for Icol=2:1:Ncol
    Columns = sprintf('%s,%s',Columns,ColCell{Icol});
end
Columns = sprintf('[%s]',Columns);


% URL = sprintf('%s%s/%s?ra=%10.6f&dec=%10.6f&radius=%10.8f&pagesize=%s&format=%s&columns=%s',...
%     InPar.BaseURL,InPar.DR,InPar.Catalog,...
%     RA,Dec,InPar.Radius,InPar.pagesize,InPar.Format,Columns);
URL = sprintf('%s%s/%s.%s?ra=%10.6f&dec=%10.6f&radius=%10.8f&pagesize=%s&columns=%s&nDetections.gte=%d',...
     InPar.BaseURL,InPar.DR,InPar.Catalog,InPar.Format,...
     RA,Dec,InPar.Radius,InPar.pagesize,Columns,InPar.nDetections);

Options = weboptions('Timeout',InPar.TimeOut);
%Data = webread('https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr1/mean?ra=210.802429&dec=54.348750&radius=0.00833333&nDetections.gte=1&pagesize=5001&format=csv');
Data = webread(URL,Options);

switch lower(InPar.OutType)
    case 'json'
        Out = Data;
    case 'mat'
        Out = Data.data;
    case 'astcat'
        Out         = AstCat;
        Out.Cat     = Data.data;
        Out.ColCell = ColCell;
        Out         = colcell2col(Out);
    otherwise
        error('Unknown OutType option');
end