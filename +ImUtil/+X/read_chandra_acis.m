function [EvtT,MetaData,Image]=read_chandra_acis(EvtFile,varargin)
% SHORT DESCRIPTION HERE
% Package: ImUtil
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: EvtT=ImUtil.X.read_chandra_acis('acis.evt2');
% Reliable: 
%--------------------------------------------------------------------------

CatField    = AstCat.CatField;
ColField    = AstCat.ColField;


DefV.ExpMap               = [];
DefV.EnergyRange          = [200 10000];
DefV.SelectGoodTimes      = true;
DefV.NinBin               = 100;

DefV.TimeColName          = 'time';
DefV.XColName             = 'x';
DefV.YColName             = 'y';
DefV.EnergyColName        = 'energy';

DefV.ImPixSize            = 0.5;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

EvtT = FITS.read_table(EvtFile);


MetaData = [];
if (InPar.SelectGoodTimes)
    [GoodInd,BadTimes,ExpTime]=ImUtil.X.find_badtimes(EvtT.(CatField),'TimeCol',EvtT.(ColField).(InPar.TimeColName),...
                                                                      'NinBin',InPar.NinBin);
    MetaData.BadTimes = BadTimes;
    MetaData.ExpTime  = ExpTime;
    EvtT.(CatField) = EvtT.(CatField)(GoodInd,:);
end

FlagE = EvtT.(CatField)(:,EvtT.(ColField).(InPar.EnergyColName))>InPar.EnergyRange(1) & ...
        EvtT.(CatField)(:,EvtT.(ColField).(InPar.EnergyColName))<InPar.EnergyRange(2);
    
EvtT.(CatField) = EvtT.(CatField)(FlagE,:);

X = EvtT.(CatField)(:,EvtT.(ColField).(InPar.XColName));
Y = EvtT.(CatField)(:,EvtT.(ColField).(InPar.YColName));

if (nargout>2)
    VecX = (min(X):InPar.ImPixSize:max(X))';
    VecY = (min(Y):InPar.ImPixSize:max(Y))';

    Image = histcounts2(X,Y,VecX,VecY);
end

