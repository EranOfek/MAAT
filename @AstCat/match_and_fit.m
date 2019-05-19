function [Res]=match_and_fit(AstC,AstR,varargin)
% Match sources in AstCat object, fit motion, and calculate rms.
% Package: @AstCat
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

DefV.natchPar             = {};
DefV.Time                 = [];
DefV.ColX                 = 'XWIN_IMAGE';
DefV.ColY                 = 'YWIN_IMAAGE';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Ncat = numel(AstC);
if isempty(InPar.Time)
    InPar.Time = (1:1:Ncat)';
end
InPar.Time = InPar.Time(:);  % make col vector

% match sources
[AstOut] = match(AstC,AstR,InPar.matchPar{:});

% prep array of matched sources
[ResM,Summary]  = astcat2matched_array(AstOut,{InPar.ColX, InPar.ColY});

Nsrc = size(ResM.(InPar.ColX),2);
Nep  = numel(InPar.Time);

% design matrix for linear motion fit
H    = [ones(Nep,1), InPar.Time];

ParX   = H \ ResM.(InPar.ColX);
ResidX = ResM.(InPar.ColX) - H*ParX;
ParY   = H \ ResM.(InPar.ColY);
ResidY = ResM.(InPar.ColY) - H*PaY;

Res.MeanPosX = nanmean(ResM.(InPar.ColX),1);
Res.MeanPosY = nanmean(ResM.(InPar.ColY),1);
Res.StdPosX  = nanstd(ResM.(InPar.ColX),1);
Res.StdPosY  = nanstd(ResM.(InPar.ColY),1);
Res.NobsX    = sum(  ~isnan(ResM.(InPar.ColX)),1 );
Res.NobsY    = sum(  ~isnan(ResM.(InPar.ColY)),1 );

Res.FitLinPosT0_X = ParX(1,:);
Res.FitLinSpeed_X = ParX(2,:);
Res.FitStd_X      = nanstd(ResidX);
Res.FitLinPosT0_Y = ParY(1,:);
Res.FitLinSpeed_Y = ParY(2,:);
Res.FitStd_Y      = nanstd(ResidY);
