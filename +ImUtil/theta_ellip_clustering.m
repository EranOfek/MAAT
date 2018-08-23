function theta_ellip_clustering(Theta,Elon,varargin)
% SHORT DESCRIPTION HERE
% Package: ImUtil
% Description: Given the orientation (THETA) and ellipticity
%              (ELONGATION; A/B) of sources, find clusters in the THETA,
%              A/B plan, and flag anomalous sources (e.g., due to aircarft
%              streaks).
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ImUtil.theta_ellip_clustering(A2.Cat(:,8),A2.Cat(:,9))
% Reliable: 
%--------------------------------------------------------------------------



DefV.MaxElon             = 3;
DefV.StepElon            = 0.2;
DefV.RangeTheta          = [-45 45];
DefV.StepTheta           = 10;
DefV.MinSN               = 10;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

ThetaEdge = (InPar.RangeTheta(1):InPar.StepTheta:InPar.RangeTheta(2)).';
ElonEdge  = (1:InPar.StepElon:InPar.MaxElon).';
N = histcounts2(Theta,Elon,ThetaEdge,ElonEdge);

MeanNumberElon = sum(N)./size(N,1);

N./MeanNumberElon

MapSN = (N./sqrt(MeanNumberElon))>InPar.MinSN

Ind   = find(sum(MapSN)>0)
Ind