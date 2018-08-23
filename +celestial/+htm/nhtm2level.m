function [Level,Nh]=nhtm2level(Nhtm)
% Given number of HTM elements calculate number of levels.
% Package: celestial
% Description: Given number of HTM elements calculate number of HTM levels.
% Input  : - Number of HTM elements.
% Output : - Number of HTM levels.
%          - Number of HTMs in the lowest level.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Level=celestial.htm.nhtm2level(43688)
% Reliable: 2
%--------------------------------------------------------------------------

Level = floor(log(Nhtm./2)./log(4));
Nh    = 2.*4.^Level;
