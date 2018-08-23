function [T,LC,Info]=get_batse_lc(TriggerNumber);
%-------------------------------------------------------------------------
% get_batse_lc function                                         Catalogue
% Description: Get a 4 chanel BATSE light curve in 64ms bins
%              (only PREB sample) from a local catalog.
% Input  : - BATSE Trigger number.
% Output : - Column vector of times in second since trigger.
%          - Four chanel light curve (chanel per column).
%            The chanels are:
%            25-55 keV, 55-110 keV, 110-320 keV, and >320 keV.
%          - Info matrix:
%            [unique BATSE trigger number,
%             total number of samples to follow  per energy channel,
%             total number of DISCLA 64-ms samples concatenated prior to PREB,
%             first PREB sample number after last 1.024-s DISCLA sample].
% Tested : Matlab 5.3
%     By : Eran O. Ofek                      Feb. 2004
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See Also: batse.dat
% Reliable: 2
%-------------------------------------------------------------------------
SubDir   = 'BATSE_LC';
BaseName = 'c64ms';
TimeStep = 0.064;  % sec
FileSep  = filesep;

FileName = sprintf('%s.%05d',BaseName,TriggerNumber);

LC       = load(FileName);
Info     = LC(1,:);
LC       = LC(2:end,:);
LC       = LC(Info(3)+1:end,:);
N        = size(LC,1);
T        = [-32.*TimeStep:TimeStep:(N-33).*TimeStep].';

