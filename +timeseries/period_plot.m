function []=period_plot(Data,Freq,varargin)
% SHORT DESCRIPTION HERE
% Package: timeseries
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV.Type                = 'Scargle';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


[PS,Peaks] = timeseries.period(Data,Freq,'Type',InPar.Type);



